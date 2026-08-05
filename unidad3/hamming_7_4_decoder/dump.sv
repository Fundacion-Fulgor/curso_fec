module dump;

initial begin
    $dumpfile("hamming_decoder.vcd");
    $dumpvars(0, hamming_decoder);
end

endmodule