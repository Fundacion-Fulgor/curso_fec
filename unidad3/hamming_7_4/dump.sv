module dump;

initial begin
    $dumpfile("hamming_7_4.vcd");
    $dumpvars(0, hamming_7_4);
end

endmodule