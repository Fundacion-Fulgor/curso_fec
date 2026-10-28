// General multiplier in GF(2^M): shift-and-add with reduction by P(x).
// Combinational; it synthesizes to an AND-XOR network.

module kes_gf_mul #(
    parameter M    = 8,
    parameter POLY = 'h1D            // P(x) without the x^M term
) (
    input  wire [M-1:0] i_a,
    input  wire [M-1:0] i_b,
    output reg  [M-1:0] o_p
);

localparam [M-1:0] P = POLY[M-1:0];

reg [M-1:0] a;                   // i_a * x^i mod P(x)
integer     i;

always @(*) begin
    o_p = {M{1'b0}};
    a   = i_a;
    for (i = 0; i < M; i = i + 1) begin
        if (i_b[i]) begin
            o_p = o_p ^ a;
        end
        a = a[M-1] ? ((a << 1) ^ P) : (a << 1);
    end
end

endmodule
