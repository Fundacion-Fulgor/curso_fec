module dump;

initial begin
    $dumpfile("hamming_encoder.vcd");
    $dumpvars(0, hamming_encoder);
end

endmodule