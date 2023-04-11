// mipstest.sv
`timescale 1ns/1ps

// Testbench for MIPS processor
module testbench();

  logic        clk, memclk;
  logic        clrn;
  
  logic [2:0] q;
  logic [31:0] a, b, alu, adr, tom, fromm, pc, ir;
  // instantiate device to be tested
  mccomp dut(clk,clrn,q,a,b,alu,adr,tom,fromm,pc,ir,memclk);
  
  // initialize test
  initial
    begin
      clrn <= 0; # 10; clrn <= 1;
    end

  // generate clock to sequence tests
  always
    begin
      clk <= 0; # 5; clk <= 1; # 5;
    end
   
  always
    begin
        memclk <= 0; # 5; memclk <= 1; # 5;     
    end


endmodule
