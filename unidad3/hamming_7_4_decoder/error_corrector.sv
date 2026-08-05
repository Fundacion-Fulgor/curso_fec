module error_corrector (
    output logic [4-1:0] o_info,
    output logic         o_is_corrected,
    input  logic         i_zero_syndrome,
    input  logic [3-1:0] i_syndrome,
    input  logic [7-1:0] i_message,
    input  logic         i_clk,
    input  logic         i_rst
);

/*
    | 1 0 1 1 1 0 0 |
H = | 1 1 1 0 0 1 0 |
    | 0 1 1 1 0 0 1 |

h0 = 110 -> 3'b011
h1 = 011 -> 3'b110
h2 = 111 -> 3'b111
h3 = 101 -> 3'b101
h4 = 100 -> 3'b001
h5 = 010 -> 3'b010
h6 = 001 -> 3'b100
*/

logic [3-1:0] error_position;
logic [7-1:0] estimated_error;
logic [7-1:0] corrected_message;
logic [4-1:0] info_d;
logic         is_corrected_d;

always_comb begin
    case (i_syndrome)
        3'b011:  error_position = 3'd0; //h0
        3'b110:  error_position = 3'd1; //h1
        3'b111:  error_position = 3'd2; //h2
        3'b101:  error_position = 3'd3; //h3
        3'b001:  error_position = 3'd4; //h4
        3'b010:  error_position = 3'd5; //h5
        3'b100:  error_position = 3'd6; //h6
        default: error_position = 3'd0;
    endcase
end

always_comb begin
    estimated_error                 = 7'd0;
    estimated_error[error_position] = 1'b1;
end

assign corrected_message = (i_zero_syndrome) ? i_message : (i_message ^ estimated_error);

always_ff @(posedge i_clk) begin
    if(i_rst) begin
        info_d         <= 4'd0;
        is_corrected_d <= 1'b0;
    end
    else begin
        info_d         <= corrected_message[3:0];
        is_corrected_d <= (i_zero_syndrome == 0);
    end
end

assign o_info         = info_d;
assign o_is_corrected = is_corrected_d;

endmodule