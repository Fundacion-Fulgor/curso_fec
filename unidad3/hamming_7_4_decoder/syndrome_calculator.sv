module syndrome_calculator (
    output logic [3-1:0] o_syndrome,
    output logic         o_zero_syndrome,
    input  logic [7-1:0] i_message,
    input  logic         i_clk,
    input  logic         i_rst
);

/*
    | 1 0 1 1 1 0 0 |
H = | 1 1 1 0 0 1 0 |
    | 0 1 1 1 0 0 1 |

      | 1 1 0 |
      | 0 1 1 |
      | 1 1 1 |
H^T = | 1 0 1 |
      | 1 0 0 |
      | 0 1 0 |
      | 0 0 1 |

s (syndrome) = r (received message) * H^T (transposed H)
r = v + e

s0 = r0 + r2 + r3 + r4
s1 = r0 + r1 + r2 + r5
s2 = r1 + r2 + r3 + r6
*/

logic [7-1:0] message_d;
logic [3-1:0] syndrome;

always_ff @(posedge i_clk) begin
    if(i_rst) begin
        message_d <= 7'd0;
    end
    else begin
        message_d <= i_message;
    end
end

always_comb begin
    syndrome[0] = message_d[0] ^ message_d[2] ^ message_d[3] ^ message_d[4];
    syndrome[1] = message_d[0] ^ message_d[1] ^ message_d[2] ^ message_d[5];
    syndrome[2] = message_d[1] ^ message_d[2] ^ message_d[3] ^ message_d[6];
end

assign o_syndrome      = syndrome;
assign o_zero_syndrome = (|syndrome == 1'b0);

// |syndrome = syndrome[0] | syndrome[1] | syndrome[2]

endmodule