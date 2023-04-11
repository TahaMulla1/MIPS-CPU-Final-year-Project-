// mipstest.sv
`timescale 1ns/1ps

// Testbench for MIPS processor
module testbench();

  logic        clk;
  logic        clrn;

  logic [31:0] pc, inst, ealu, malu, wdi, mb;
  logic mwmem;

  // instantiate device to be tested
  pipelinedcpu dut(clk,clrn,pc,inst,ealu,malu,mb,wdi,mwmem);
  
  // initialize test
  initial
    begin
      clrn <= 0; # 50; clrn <= 1;
    end

  // generate clock to sequence tests
  always
    begin
      clk <= 1; # 5; clk <= 0; # 5;
    end

  // check results
  always @(negedge clk)
    begin
       if(pc === 32'h0000008c & mb === 32'd2) begin
         $display("Simulation succeeded");
         @(negedge clk);
         $stop;
        end else if (pc === 32'h0000008c & mb != 32'd2) begin
          $display("Simulation failed");
          $stop;
      end
    end
endmodule
