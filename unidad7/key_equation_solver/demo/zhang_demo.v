// key_equation_solver on Example 31 of Zhang, "VLSI Architectures for Modern Error-Correcting
// Codes": RS(7, 3) over GF(2^3), P(x) = x^3 + x + 1, t = 2, S = (0, a, a^2, a^6).
// Every coefficient of the registers gets its own 3-bit signal, so the waveform viewer can show
// it as a power of alpha with the translate filter alpha_gf8.txt.

module zhang_demo (
    input  wire        i_clk,
    input  wire        i_rst,
    input  wire        i_valid,
    input  wire [11:0] i_syn,
    output wire        o_valid,
    output wire [8:0]  o_lambda,
    output wire [5:0]  o_omega
);

    key_equation_solver #(.M(3), .POLY('h3), .T(2)) u_kes (
        .i_clk    (i_clk),
        .i_rst    (i_rst),
        .i_valid  (i_valid),
        .i_syn    (i_syn),
        .o_valid  (o_valid),
        .o_lambda (o_lambda),
        .o_omega  (o_omega)
    );

    // control: iteration r, k, gamma and MC
    wire       busy  = u_kes.u_ctrl.busy_q;
    wire [1:0] r     = u_kes.u_ctrl.count_q;
    wire [3:0] k     = u_kes.u_ctrl.k_q;
    wire [2:0] gamma = u_kes.u_ctrl.gamma_q;
    wire       mc    = u_kes.mc;

    // rDC block: delta-hat_i and theta-hat_i
    wire [2:0] delta_0 = u_kes.delta_q[2:0];
    wire [2:0] delta_1 = u_kes.delta_q[5:3];
    wire [2:0] delta_2 = u_kes.delta_q[8:6];
    wire [2:0] delta_3 = u_kes.delta_q[11:9];
    wire [2:0] theta_0 = u_kes.theta_q[2:0];
    wire [2:0] theta_1 = u_kes.theta_q[5:3];
    wire [2:0] theta_2 = u_kes.theta_q[8:6];
    wire [2:0] theta_3 = u_kes.theta_q[11:9];

    // ELU block: lambda_i and b_i
    wire [2:0] lambda_0 = u_kes.lam_q[2:0];
    wire [2:0] lambda_1 = u_kes.lam_q[5:3];
    wire [2:0] lambda_2 = u_kes.lam_q[8:6];
    wire [2:0] b_0      = u_kes.b_q[2:0];
    wire [2:0] b_1      = u_kes.b_q[5:3];
    wire [2:0] b_2      = u_kes.b_q[8:6];

    // outputs
    wire [2:0] out_lambda_0 = o_lambda[2:0];
    wire [2:0] out_lambda_1 = o_lambda[5:3];
    wire [2:0] out_lambda_2 = o_lambda[8:6];
    wire [2:0] out_omega_0  = o_omega[2:0];
    wire [2:0] out_omega_1  = o_omega[5:3];

endmodule
