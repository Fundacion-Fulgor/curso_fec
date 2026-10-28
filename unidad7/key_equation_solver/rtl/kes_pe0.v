// PE0 of the ELU block (Sarwate and Shanbhag, Fig. 3a), without its latches:
//   lambda_i(r+1) = gamma(r) * lambda_i(r) + delta_0(r) * b_(i-1)(r)
//   b_i(r+1)      = MC(r) ? lambda_i(r) : b_(i-1)(r)

module kes_pe0 #(
    parameter M    = 8,
    parameter POLY = 'h1D
) (
    input  wire [M-1:0] i_gamma,
    input  wire [M-1:0] i_d0,        // delta_0(r), the discrepancy
    input  wire         i_mc,
    input  wire [M-1:0] i_lambda,    // lambda_i(r)
    input  wire [M-1:0] i_b,         // b_(i-1)(r), from the left neighbour
    output wire [M-1:0] o_lambda,    // lambda_i(r+1)
    output wire [M-1:0] o_b          // b_i(r+1)
);

wire [M-1:0] p_gamma, p_d0;

kes_gf_mul #(.M(M), .POLY(POLY)) u_mul_gamma (.i_a(i_gamma), .i_b(i_lambda), .o_p(p_gamma));
kes_gf_mul #(.M(M), .POLY(POLY)) u_mul_d0    (.i_a(i_d0),    .i_b(i_b),      .o_p(p_d0));

assign o_lambda = p_gamma ^ p_d0;
assign o_b      = i_mc ? i_lambda : i_b;

endmodule
