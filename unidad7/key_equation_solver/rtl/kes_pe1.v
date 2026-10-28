// PE1 of the rDC block (Sarwate and Shanbhag, Fig. 4a), without its latches:
//   delta_i(r+1) = gamma(r) * delta_(i+1)(r) + delta_0(r) * theta_i(r)
//   theta_i(r+1) = MC(r) ? delta_(i+1)(r) : theta_i(r)

module kes_pe1 #(
    parameter M    = 8,
    parameter POLY = 'h1D
) (
    input  wire [M-1:0] i_gamma,
    input  wire [M-1:0] i_d0,        // delta_0(r), the discrepancy
    input  wire         i_mc,
    input  wire [M-1:0] i_delta,     // delta_(i+1)(r), from the right neighbour
    input  wire [M-1:0] i_theta,     // theta_i(r)
    output wire [M-1:0] o_delta,     // delta_i(r+1)
    output wire [M-1:0] o_theta      // theta_i(r+1)
);

wire [M-1:0] p_gamma, p_d0;

kes_gf_mul #(.M(M), .POLY(POLY)) u_mul_gamma (.i_a(i_gamma), .i_b(i_delta), .o_p(p_gamma));
kes_gf_mul #(.M(M), .POLY(POLY)) u_mul_d0    (.i_a(i_d0),    .i_b(i_theta), .o_p(p_d0));

assign o_delta = p_gamma ^ p_d0;
assign o_theta = i_mc ? i_delta : i_theta;

endmodule
