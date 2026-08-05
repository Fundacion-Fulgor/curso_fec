module hamming_7_4 (
    output logic [4-1:0] o_info,
    output logic         o_is_corrected,
    input  logic [4-1:0] i_info,
    input  logic [7-1:0] i_error_vector,
    input  logic         i_clk,
    input  logic         i_rst
);

/*
    i_info ---------> [reg] ---> encoder ---> (+) ---> decoder ---> o_info
                                               ^
    i_error_vector -> [reg] --------------------

    el vector de error se registra junto con la informacion, asi la suma
    queda alineada en latencia
*/

logic [4-1:0] info_d;
logic [7-1:0] error_vector_d;
logic [7-1:0] message;
logic [7-1:0] received;

always_ff @(posedge i_clk) begin
    if(i_rst) begin
        info_d         <= 'd0;
        error_vector_d <= 'd0;
    end
    else begin
        info_d         <= i_info;
        error_vector_d <= i_error_vector;
    end
end

hamming_encoder u_hamming_encoder(
    .o_message(message),
    .i_info(info_d)
);

assign received = message ^ error_vector_d;

hamming_decoder u_hamming_decoder(
    .o_info(o_info),
    .o_is_corrected(o_is_corrected),
    .i_message(received),
    .i_clk(i_clk),
    .i_rst(i_rst)
);

endmodule