// Key equation solver for a t-error-correcting Reed-Solomon code over GF(2^M):
// reformulated inversionless Berlekamp-Massey (riBM) of Sarwate and Shanbhag,
// "High-speed architectures for Reed-Solomon decoders", IEEE TVLSI, 2001, Fig. 5.
//
//   rDC block: 2T PE1 holding delta_i and theta_i (i = 0 .. 2T-1)
//   ELU block: T+1 PE0 holding lambda_i and b_i  (i = 0 .. T)
//
// Buses hold the coefficient of degree i in bits [M*i+M-1 : M*i].
//   i_syn    S_0 .. S_(2T-1), with S_j = r(alpha^j)
//   o_lambda Lambda(x), scaled by a non-zero constant (lambda_0 is not 1)
//   o_omega  Omega_h(x): the coefficients of x^(2T) .. x^(3T-1) of Lambda(x) S(x)
// Timing: docs/kes_protocol.md. Latency and interval are 2T clocks.

module key_equation_solver #(
    parameter M    = 8,
    parameter POLY = 'h1D,
    parameter T    = 8
) (
    input  wire               i_clk,
    input  wire               i_rst,
    input  wire               i_valid,
    input  wire [2*T*M-1:0]   i_syn,
    output reg                o_valid,
    output reg  [(T+1)*M-1:0] o_lambda,
    output reg  [T*M-1:0]     o_omega
);

localparam DW = 2 * T * M;                   // rDC width
localparam LW = (T + 1) * M;                 // ELU width

// registers, named as the variables of the riBM pseudocode
reg  [DW-1:0] delta_q, theta_q;
reg  [LW-1:0] lam_q, b_q;

// next state
wire [DW-1:0] delta_d, theta_d;
wire [LW-1:0] lam_d, b_d;

wire          step, last, mc;
wire [M-1:0]  gamma;
wire [M-1:0]  d0 = delta_q[M-1:0];           // the discrepancy is always delta_0

wire [DW-1:0] delta_right = {{M{1'b0}}, delta_q[DW-1:M]};   // delta_(i+1), delta_2T = 0
wire [LW-1:0] b_left      = {b_q[LW-M-1:0], {M{1'b0}}};     // b_(i-1),     b_-1 = 0

kes_ctrl #(.M(M), .T(T)) u_ctrl (
    .i_clk   (i_clk),
    .i_rst   (i_rst),
    .i_valid (i_valid),
    .i_d0    (d0),
    .o_step  (step),
    .o_last  (last),
    .o_mc    (mc),
    .o_gamma (gamma)
);

genvar i;
generate
    for (i = 0; i < 2 * T; i = i + 1) begin : g_rdc
        kes_pe1 #(.M(M), .POLY(POLY)) u_pe1 (
            .i_gamma (gamma),
            .i_d0    (d0),
            .i_mc    (mc),
            .i_delta (delta_right[M*i +: M]),
            .i_theta (theta_q[M*i +: M]),
            .o_delta (delta_d[M*i +: M]),
            .o_theta (theta_d[M*i +: M])
        );
    end
    for (i = 0; i <= T; i = i + 1) begin : g_elu
        kes_pe0 #(.M(M), .POLY(POLY)) u_pe0 (
            .i_gamma  (gamma),
            .i_d0     (d0),
            .i_mc     (mc),
            .i_lambda (lam_q[M*i +: M]),
            .i_b      (b_left[M*i +: M]),
            .o_lambda (lam_d[M*i +: M]),
            .o_b      (b_d[M*i +: M])
        );
    end
endgenerate

// datapath: load the syndromes, iterate while busy, hold when idle
always @(posedge i_clk) begin
    if (i_valid) begin
        delta_q <= i_syn;
        theta_q <= i_syn;
        lam_q   <= {{(LW-1){1'b0}}, 1'b1};
        b_q     <= {{(LW-1){1'b0}}, 1'b1};
    end
    else if (step) begin
        delta_q <= delta_d;
        theta_q <= theta_d;
        lam_q   <= lam_d;
        b_q     <= b_d;
    end
end

// outputs: the result of the last iteration
always @(posedge i_clk) begin
    if (i_rst) begin
        o_valid  <= 1'b0;
        o_lambda <= {LW{1'b0}};
        o_omega  <= {(T*M){1'b0}};
    end
    else begin
        o_valid <= last;
        if (last) begin
            o_lambda <= lam_d;
            o_omega  <= delta_d[T*M-1:0];
        end
    end
end

endmodule
