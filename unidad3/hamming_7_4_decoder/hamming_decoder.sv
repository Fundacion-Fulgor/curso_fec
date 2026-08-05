module hamming_decoder (
    output logic [4-1:0] o_info,
    output logic         o_is_corrected,
    input  logic [7-1:0] i_message,
    input  logic         i_clk,
    input  logic         i_rst
);

logic [7-1:0] message_d;
logic [3-1:0] syndrome;
logic         zero_syndrome;

always_ff @(posedge i_clk) begin
    message_d <= i_message;
end

syndrome_calculator u_syndrome_calculator(
    .o_syndrome      (syndrome),
    .o_zero_syndrome (zero_syndrome),
    .i_message       (i_message),
    .i_clk           (i_clk),
    .i_rst           (i_rst)
);

error_corrector u_error_corrector(
    .o_info          (o_info),
    .o_is_corrected  (o_is_corrected),
    .i_zero_syndrome (zero_syndrome),
    .i_syndrome      (syndrome),
    .i_message       (message_d),
    .i_clk           (i_clk),
    .i_rst           (i_rst)
);
    
endmodule