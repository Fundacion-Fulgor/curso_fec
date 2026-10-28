// Smoke test DUT: registered XOR of two bytes.
module smoke_dut (
    input  wire       i_clk,
    input  wire       i_rst,
    input  wire [7:0] i_a,
    input  wire [7:0] i_b,
    output reg  [7:0] o_q
);
    always @(posedge i_clk) begin
        if (i_rst) o_q <= 8'h00;
        else       o_q <= i_a ^ i_b;
    end
endmodule
