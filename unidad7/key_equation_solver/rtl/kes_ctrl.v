// Control unit (Sarwate and Shanbhag, Fig. 2) plus the sequencing of the protocol:
// a word is loaded with i_valid, takes one iteration per clock while busy, and the
// iteration with count = 2T - 1 is the last one. Loading always wins and restarts.
//
//   MC(r)   = delta_0(r) != 0 and k(r) >= 0
//   gamma   <- MC ? delta_0 : gamma
//   k       <- MC ? -k - 1  : k + 1        (-k - 1 is ~k in two's complement)

module kes_ctrl #(
    parameter M = 8,
    parameter T = 8
) (
    input  wire         i_clk,
    input  wire         i_rst,
    input  wire         i_valid,
    input  wire [M-1:0] i_d0,        // delta_0(r)
    output wire         o_step,      // an iteration happens on this edge
    output wire         o_last,      // and it is the last one of the word
    output wire         o_mc,
    output wire [M-1:0] o_gamma
);

localparam N  = 2 * T;                       // iterations per word
localparam CW = (N > 2) ? $clog2(N) : 1;     // iteration counter width
localparam KW = $clog2(2 * T + 1) + 1;       // k in [-2T, 2T], signed

localparam          NM1  = N - 1;
localparam [CW-1:0] LAST = NM1[CW-1:0];

reg          busy_q;
reg [CW-1:0] count_q;                        // r, the current state
reg [M-1:0]  gamma_q;
reg [KW-1:0] k_q;

assign o_step  = busy_q;
assign o_last  = busy_q && (count_q == LAST);
assign o_mc    = (|i_d0) && !k_q[KW-1];
assign o_gamma = gamma_q;

always @(posedge i_clk) begin
    if (i_rst) begin
        busy_q  <= 1'b0;
        count_q <= {CW{1'b0}};
    end
    else if (i_valid) begin
        busy_q  <= 1'b1;
        count_q <= {CW{1'b0}};
    end
    else if (busy_q) begin
        busy_q  <= !o_last;
        count_q <= count_q + 1'b1;
    end
end

always @(posedge i_clk) begin
    if (i_valid) begin
        gamma_q <= {{(M-1){1'b0}}, 1'b1};
        k_q     <= {KW{1'b0}};
    end
    else if (busy_q) begin
        if (o_mc) begin
            gamma_q <= i_d0;
            k_q     <= ~k_q;
        end
        else begin
            k_q     <= k_q + 1'b1;
        end
    end
end

endmodule
