// mipssingle.sv

// Single-cycle MIPS processor

module mips(input  logic        clk, reset,
            output logic [31:0] pc,
            input  logic [31:0] instr,
            output logic        memwrite,
            output logic [31:0] aluout, writedata,
            input  logic [31:0] readdata);

  logic       memtoreg, alusrc, regdst, 
              regwrite, jump, pcsrc, zero, jump_reg, zero_imm, jal;  //adding two new control signals, jump_reg and zero_imm
  logic [3:0] alucontrol;                                       //jump_reg for jr and zero_imm for zero extend
 //alucontrol changed to 4 bits to accomodate all the functions
 //also passing the two new control signals into controller and datapath
  controller c(instr[31:26], instr[5:0], zero,
               memtoreg, memwrite, pcsrc,
               alusrc, regdst, regwrite, jump,
               jump_reg, zero_imm, jal, alucontrol);

  datapath dp(instr[31:26], clk, reset, memtoreg, pcsrc,
              alusrc, regdst, regwrite, jump,
              jump_reg, zero_imm, jal, alucontrol,
              zero, pc, instr,
              aluout, writedata, readdata);
endmodule

module controller(input  logic [5:0] op, funct,
                  input  logic       zero,
                  output logic       memtoreg, memwrite,
                  output logic       pcsrc, alusrc,
                  output logic       regdst, regwrite,
                  output logic       jump, jump_reg, zero_imm, jal,
                  output logic [3:0] alucontrol);

  logic [2:0] aluop;        //aluop is changed to 3 bits
  logic       branch;
//also funct and the two new control signals are passed into main decoder,
//function is passed for JR since it requires a different set of control signals configuration,
//despite being R-type
  maindec md(op, funct, memtoreg, memwrite, branch,
             alusrc, regdst, regwrite, jump, jump_reg, zero_imm, jal, aluop);
  aludec  ad(funct, aluop, alucontrol);

  assign pcsrc = branch & zero;
endmodule

module maindec(input  logic [5:0] op, funct,
               output logic       memtoreg, memwrite,
               output logic       branch, alusrc,
               output logic       regdst, regwrite,
               output logic       jump, jump_reg, zero_imm, jal,
               output logic [2:0] aluop);

  logic [12:0] controls;    //Control logic changed to 12 bits because of the increase in aluop
                            //and the two new control signals
  assign {regwrite, regdst, alusrc, branch, memwrite,
          memtoreg, jump, jump_reg, zero_imm, jal, aluop} = controls;

  always_comb
    case(op)
    //for the R-type instruction, since jr is an exception and requires a different set of signals,
    //i use function to differentiate between them and pass the right control signals into JR
      6'b000000: begin
        if(funct == 6'b001000)
          controls <= 13'b00000001x0000;        //JR (jump_reg is 1)
        else
          controls <= 13'b11000000x0010; // RTYPE
      end
      6'b100011: controls <= 13'b1010010000000; // LW
      6'b101011: controls <= 13'b0010100000000; // SW
      6'b000100: controls <= 13'b00010000x0001; // BEQ
      6'b001000: controls <= 13'b1010000000000; // ADDI
      6'b000010: controls <= 13'b00000010x0000; // J
      6'b000011: controls <= 13'b00000010x1000; // Jal
      6'b001110: controls <= 13'b1010000010011; //xori (zero_imm is 1)
      6'b001100: controls <= 13'b1010000010100; //andi (zero_imm is 1)
      6'b001101: controls <= 13'b1010000010101; //ori  (zero_imm is 1)
      6'b000111: controls <= 13'b00010000x0000; // BGTZ (branch is 1)
      default:   controls <= 13'bxxxxxxxxxxxxx; // illegal op
    endcase
endmodule

module aludec(input  logic [5:0] funct,
              input  logic [2:0] aluop,
              output logic [3:0] alucontrol);
//as per the theory, the alucontrol is assigned depending on the aluop or function
  always_comb
    case(aluop)
      3'b000: alucontrol <= 4'b0010;  // add (for lw/sw/addi)
      3'b001: alucontrol <= 4'b1010;  // sub (for beq)
      3'b011: alucontrol <= 4'b0101;  // xor (for xori)
      3'b100: alucontrol <= 4'b0000;  // and (for andi)
      3'b101: alucontrol <= 4'b0001;  // or (ori)
      default: case(funct)          // R-type instructions
          6'b100010: alucontrol <= 4'b1010; // sub
          6'b100100: alucontrol <= 4'b0000; // and
          6'b100000: alucontrol <= 4'b0010; // add
          6'b100101: alucontrol <= 4'b0001; // or
          6'b101010: alucontrol <= 4'b1011; // slt
          6'b000010: alucontrol <= 4'b0100; // srl
          6'b100110: alucontrol <= 4'b0101; // xor
          6'b000110: alucontrol <= 4'b0110; // srlv
          default:   alucontrol <= 4'bxxxx; // ???
        endcase
    endcase      
endmodule

module datapath(input  logic [5:0]  op,
                input  logic        clk, reset,
                input  logic        memtoreg, pcsrc,
                input  logic        alusrc, regdst,
                input  logic        regwrite, jump, jump_reg, zero_imm, jal, //the new control logic 
                input  logic [3:0]  alucontrol,
                output logic        zero,
                output logic [31:0] pc,
                input  logic [31:0] instr,
                output logic [31:0] aluout, writedata,
                input  logic [31:0] readdata);

  logic [4:0]  beforejalreg, writereg;
  logic [31:0] pcnext, pcnext_reg, pcnextbr, pcplus4, pcbranch; //create a new logic pcnext_reg for JR
  logic [31:0] signimm, signimmsh, zeimm, immediate;        //create zeimm and immediate for the choice between zero_imm and sign_imm
  logic [31:0] srca, srcb;
  logic [31:0] before_jal_result, result;
  logic branch_test;   // 1-bit logic to determine whether bgtz or beq is selected
 //Combinational block to determine if the instruction is bgtz or beq
  always_comb begin
    if(op == 6'b000111)             //If op == 000111 then it is bgtz and if(a>0) is checked
      branch_test = (!(srca[31] == 1'b1)) && (srca > 32'b0);
    else
      branch_test = pcsrc;           //or else it is branch equal
  end

  // next PC logic
  flopr #(32) pcreg(clk, reset, pcnext, pc);
  adder       pcadd1(pc, 32'b100, pcplus4);
  sl2         immsh(signimm, signimmsh);
  adder       pcadd2(pcplus4, signimmsh, pcbranch);
  mux2 #(32)  pcbrmux(pcplus4, pcbranch, branch_test, pcnextbr);    //pcsrc replaced by branch_test since now
  mux2 #(32)  pcmux(pcnextbr, {pcplus4[31:28],                      //that is the selection signal
                    instr[25:0], 2'b00}, jump, pcnext_reg);
  mux2 #(32)  pcnextreg(pcnext_reg, srca, jump_reg, pcnext);        //MUX for jr, whether rs gets passed into PC or not

  // register file logic
  regfile     rf(clk, regwrite, instr[25:21], instr[20:16], 
                 writereg, result, srca, writedata);
  mux2 #(5)   wrmux(instr[20:16], instr[15:11],
                    regdst, beforejalreg);
  mux2 #(5)   rjalmux(beforejalreg, 5'd31, jal, writereg);
  mux2 #(32)  resmux(aluout, readdata, memtoreg, before_jal_result);
  mux2 #(32)  thirtyonemux(before_jal_result, pcplus4, jal, result);
  signext     se(instr[15:0], signimm);
  zeroimm     ze(instr[15:0], zeimm);       // to make a zero immediate value
  mux2 #(32)  imm_mux(signimm, zeimm, zero_imm, immediate); // to choose between a zero immediate value and sign immediate value
  // ALU logic
  mux2 #(32)  srcbmux(writedata, immediate, alusrc, srcb);  //immediate is the result of the above choice
  alu         alu(srca, srcb, instr[10:6], alucontrol, aluout, zero); //instr[10:6] passed into ALU since that contains the shamt
endmodule

module regfile(input  logic        clk, 
               input  logic        we3, 
               input  logic [4:0]  ra1, ra2, wa3, 
               input  logic [31:0] wd3, 
               output logic [31:0] rd1, rd2);

  logic [31:0] rf[31:0];

  // three ported register file
  // read two ports combinationally
  // write third port on rising edge of clk
  // register 0 hardwired to 0
  // note: for pipelined processor, write third port
  // on falling edge of clk

  always_ff @(posedge clk)
    if (we3) rf[wa3] <= wd3;	

  assign rd1 = (ra1 != 0) ? rf[ra1] : 0;
  assign rd2 = (ra2 != 0) ? rf[ra2] : 0;
endmodule

module adder(input  logic [31:0] a, b,
             output logic [31:0] y);

  assign y = a + b;
endmodule

module sl2(input  logic [31:0] a,
           output logic [31:0] y);

  // shift left by 2
  assign y = {a[29:0], 2'b00};
endmodule

module signext(input  logic [15:0] a,
               output logic [31:0] y);
              
  assign y = {{16{a[15]}}, a};
endmodule

module zeroimm(input  logic [15:0] a,
               output logic [31:0] y);
              
  assign y = {16'b0, a};
endmodule

module flopr #(parameter WIDTH = 8)
              (input  logic             clk, reset,
               input  logic [WIDTH-1:0] d, 
               output logic [WIDTH-1:0] q);

  always_ff @(posedge clk, posedge reset)
    if (reset) q <= 0;
    else       q <= d;
endmodule

module mux2 #(parameter WIDTH = 8)
             (input  logic [WIDTH-1:0] d0, d1, 
              input  logic             s, 
              output logic [WIDTH-1:0] y);

  assign y = s ? d1 : d0; 
endmodule

module alu(input  logic [31:0] a, b,
           input  logic [4:0]  shams,
           input  logic [3:0]  alucontrol,
           output logic [31:0] result,
           output logic        zero);

  logic [31:0] condinvb, sum;

  assign condinvb = alucontrol[3] ? ~b : b;
  addsub32 as32 (sum, a, condinvb, alucontrol[3]);
  //assign sum = a + condinvb + alucontrol[3];

  always_comb
    case (alucontrol[2:0])
      3'b000: result = a & condinvb;        //AND
      3'b001: result = a | condinvb;        //OR
      3'b010: result = sum;         //SUM
      3'b011: result = sum[31];     //SLT
      3'b100: result = condinvb >> shams;       //SRL
      3'b101: result = a ^ condinvb;        //XOR
      3'b110: result = condinvb >> a[4:0];  //SRLV
    endcase

  assign zero = (result == 32'b0);
endmodule
