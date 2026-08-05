module hamming_encoder (
    output logic [7-1:0] o_message,
    input  logic [4-1:0] i_info
);

/*
    | 1 0 0 0 1 1 0 |
G = | 0 1 0 0 0 1 1 |
    | 0 0 1 0 1 1 1 |
    | 0 0 0 1 1 0 1 |

v (message) = u (info) * G (generator matrix)
*/

logic [3-1:0] parity;
logic [7-1:0] message;

assign message[0] = i_info[0];
assign message[1] = i_info[1];
assign message[2] = i_info[2];
assign message[3] = i_info[3];

assign parity[0] = i_info[0] ^ i_info[2] ^ i_info[3];
assign parity[1] = i_info[0] ^ i_info[1] ^ i_info[2];
assign parity[2] = i_info[1] ^ i_info[2] ^ i_info[3];

assign message[4] = parity[0];
assign message[5] = parity[1];
assign message[6] = parity[2];

assign o_message = message;

endmodule