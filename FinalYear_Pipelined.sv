module pipelinedcpu (input logic clk,clrn, output logic [31:0] pc,inst,ealu,malu,mb,wdi, output logic mwmem); // pipelined cpu
  /*
  clk, clrn; // clock and reset // plus inst mem
  pc; // program counter // and data mem
  inst; // instruction in ID stage
  ealu; // alu result in EXE stage
  malu; // alu result in MEM stage
  mb;  // operand b in MEM stage
  wdi; // data to be written into register file
  mwmem; // memory write in MEM stage 
  */

  // signals in IF stage
  logic [31:0] pc4, ins, npc; // pc+4 in IF stage // instruction in IF stage // next pc in IF stage
 
  // signals in ID stage
  logic [31:0] jpc, bpc, dpc4, da,db, dimm; //jump pc // branch pc // pc+4 in ID stage // two operands a and b in ID stage // 32-bit extended immediate in ID stage
  logic [4:0]  drn; // destination register number in ID stage
  logic [4:0]  daluc; // alu control in ID stage
  logic [1:0]  pcsrc; // next pc (npc) select in ID stage
  logic wpcir, dwreg, dm2reg, dwmem, daluimm, dshift, djal; // pipepc and pipeir write enable // register file write enable in ID stage // memory to register in ID stage // memory write in ID stage
                                                            // alu input b is an immediate in ID stage // shift in ID stage // jal in ID stage

  // signals in EXE stage
  logic [31:0] epc4, ea,eb, eimm; // pc+4 in EXE stage // two operands a and b in EXE stage // 32-bit extended immediate in EXE stage
  logic [4:0] ern0, ern; // temporary register number in WB stage // destination register number in EXE stage
  logic [4:0] ealuc; // alu control in EXE stage
  logic ewreg, em2reg, ewmem, ealuimm, eshift, ejal; // register file write enable in EXE stage // memory to register in EXE stage // memory write in EXE stage // alu input b is an immediate in EXE stage
                                                     // shift in EXE stage // jal in EXE stage

  // signals in MEM stage
  logic [31:0] mmo;  // memory data out in MEM stage
  logic [4:0] mrn; // destination register number in MEM stage
  logic mwreg, mm2reg; // register file write enable in MEM stage // memory to register in MEM stage

  // signals in WB stage
  logic [31:0] wmo, walu; // memory data out in WB stage // alu result in WB stage
  logic [4:0]  wrn; // destination register number in WB stage
  logic wwreg, wm2reg; // register file write enable in WB stage // memory to register in WB stage

  // program counter
  pipepc prog_cnt (npc,wpcir,clk,clrn,pc);
  pipeif if_stage (pcsrc,pc,bpc,da,jpc,npc,pc4,ins); // IF stage
  // IF/ID pipeline register
  pipeir fd_reg (pc4,ins,wpcir,clk,clrn,dpc4,inst);
  pipeid id_stage (mwreg,mrn,ern,ewreg,em2reg,mm2reg,dpc4,inst,wrn,wdi,
                  ealu,malu,mmo,wwreg,clk,clrn,bpc,jpc,pcsrc,wpcir,
                  dwreg,dm2reg,dwmem,daluc,daluimm,da,db,dimm,drn,
                  dshift,djal); // ID stage
  // ID/EXE pipeline register
  pipedereg de_reg (dwreg,dm2reg,dwmem,daluc,daluimm,da,db,dimm,drn,
                    dshift,djal,dpc4,clk,clrn,ewreg,em2reg,ewmem,
                    ealuc,ealuimm,ea,eb,eimm,ern0,eshift,ejal,epc4);
  pipeexe exe_stage (ealuc,ealuimm,ea,eb,eimm,eshift,ern0,epc4,ejal,
                    ern,ealu); // EXE stage
  // EXE/MEM pipeline register
  pipeemreg em_reg (ewreg,em2reg,ewmem,ealu,eb,ern,clk,clrn,mwreg,
                    mm2reg,mwmem,malu,mb,mrn);
  pipemem mem_stage (mwmem,malu,mb,clk,mmo); // MEM stage
  // MEM/WB pipeline register
  pipemwreg mw_reg (mwreg,mm2reg,mmo,malu,mrn,clk,clrn,wwreg,wm2reg,
                    wmo,walu,wrn);
  pipewb wb_stage (walu,wmo,wm2reg,wdi); // WB stage
endmodule

module pipepc (input logic [31:0] npc, input logic wpc, clk, clrn, output logic [31:0] pc); // program counter // clock and reset // pc write enable // next pc // program counter
  // dffe32 (d, clk,clrn,e, q);
  dffe32 prog_cntr (npc,clk,clrn,wpc,pc); // program counter
endmodule

module pipeif (input logic [1:0] pcsrc, input logic [31:0] pc, bpc, rpc, jpc, output logic [31:0] npc, pc4, ins); // IF stage // next pc (npc) select // program counter // branch target // jump target of jr 
                                                                                                             // jump target of j/jal // next pc // pc + 4 // inst from inst mem
  mux4 #(32) next_pc (pc4,bpc,rpc,jpc,pcsrc,npc); // npc select
  cla32 pc_plus4 (pc4,pc,32'd4,1'b0); // pc + 4
  pl_inst_mem inst_mem (pc,ins); // inst mem
endmodule

module pipeir (input logic [31:0] pc4, ins, input logic wir, clk, clrn, output logic [31:0] dpc4, inst); // IF/ID pipeline register // clock and reset // write enable // pc + 4 in IF stage
                                                                                                         // instruction in IF stage // pc + 4 in ID stage // instruction in ID stage
  // dffe32 (d, clk,clrn,e, q);
  dffe32 pc_plus4 (pc4,clk,clrn,wir,dpc4); // pc+4 register
  dffe32 instruction (ins,clk,clrn,wir,inst); // inst register
endmodule

module pipeid (input logic mwreg, input logic [4:0] mrn, ern, input logic ewreg, em2reg, mm2reg, input logic [31:0] dpc4, inst, input logic [4:0] wrn, input logic [31:0] wdi, ealu,
              malu, mmo, input logic wwreg, clk, clrn, output logic [31:0] bpc, jpc, output logic [1:0] pcsrc, output logic nostall, wreg, m2reg,
              wmem, output logic [4:0] aluc, output logic aluimm, output logic [31:0] a, b, dimm, output logic [4:0] rn, output logic shift, jal);// ID stage

  logic [5:0] op, func;
  logic [4:0] rs, rt, rd;
  logic [15:0] imm, s16;
  logic [25:0] addr;
  logic [31:0] qa, qb, dis;
  logic [1:0]  fwda, fwdb;
  logic        regrt, sext, rsrtequ, bgtz_flag;
  
  //using initial
  /*
  initial
  begin
    op = inst[31:26];
    func = inst[5:0];
    rs = inst[25:21];
    rt = inst[20:16]; // rt
    rd = inst[15:11]; // rd
    imm = inst[15:0]; // immediate
    s16 = {16{sext & inst[15]}}; // 16-bit signs
    addr = inst[25:0]; // address
    dis = {dimm[29:0],2'b00}; // branch offset
    rsrtequ = ~|(a^b);      // reg[rs] == reg[rt]
  end
  */
  
  //using assign
  assign op = inst[31:26];
  assign func = inst[5:0];
  assign rs = inst[25:21];
  assign rt = inst[20:16]; // rt
  assign rd = inst[15:11]; // rd
  assign imm = inst[15:0]; // immediate
  assign s16 = {16{sext & inst[15]}}; // 16-bit signs
  assign addr = inst[25:0]; // address
  assign dis = {dimm[29:0],2'b00}; // branch offset
  assign rsrtequ = ~|(a^b);      // reg[rs] == reg[rt]
  assign bgtz_flag = (!(a[31] == 1'b1)) && (a > 32'b0); //to check if greater than zero and also not zero

  pipeidcu cu (mwreg,mrn,ern,ewreg,em2reg,mm2reg, // control unit
              rsrtequ,bgtz_flag,func,op,rs,rt,wreg,m2reg,
              wmem,aluc,regrt,aluimm,fwda,fwdb,
              nostall,sext,pcsrc,shift,jal);
  regfile r_f (rs,rt,wdi,wrn,wwreg,~clk,clrn,qa,qb); // register file
  mux2 #(5) d_r (rd,rt,regrt,rn); // select dest reg #
  mux4 #(32) s_a (qa,ealu,malu,mmo,fwda,a); // forward for a
  mux4 #(32) s_b (qb,ealu,malu,mmo,fwdb,b); // forward for b
  cla32 b_adr (bpc,dpc4,dis,1'b0); // branch target
  //adder bpcadd (dpc4,dis,bpc);
  assign dimm = {s16,imm}; // 32-bit imm
  assign jpc = {dpc4[31:28],addr,2'b00}; // jump target
endmodule

module pipeidcu (input logic mwreg, input logic [4:0] mrn, ern, input logic ewreg, em2reg, mm2reg, rsrtequ, bgtz_flag, input logic [5:0] func, op, input logic [4:0] rs, rt,
                output logic wreg, m2reg, wmem, output logic [4:0] aluc, output logic regrt, aluimm, output logic [1:0] fwda, fwdb, output logic nostall, sext,
                output logic [1:0] pcsrc, output logic shift, jal); // control unit in ID stage
  // instruction decode
  logic rtype, i_add, i_sub, i_and, i_or, i_xor, i_sll, i_srl, i_sra, i_jr, i_srlv, i_slt;
  logic i_addi, i_andi, i_ori, i_xori, i_lw, i_sw, i_beq, i_bne, i_bgtz, i_lui, i_j, i_jal;

  and (rtype,~op[5],~op[4],~op[3],~op[2],~op[1],~op[0]); // r format
  and (i_add,rtype, func[5],~func[4],~func[3],~func[2],~func[1],~func[0]);
  and (i_sub,rtype, func[5],~func[4],~func[3],~func[2], func[1],~func[0]);
  and (i_and,rtype, func[5],~func[4],~func[3], func[2],~func[1],~func[0]);
  and (i_or, rtype, func[5],~func[4],~func[3], func[2],~func[1], func[0]);
  and (i_xor,rtype, func[5],~func[4],~func[3], func[2], func[1],~func[0]);
  and (i_sll,rtype,~func[5],~func[4],~func[3],~func[2],~func[1],~func[0]);
  and (i_srl,rtype,~func[5],~func[4],~func[3],~func[2], func[1],~func[0]);
  and (i_sra,rtype,~func[5],~func[4],~func[3],~func[2], func[1], func[0]);
  and (i_jr, rtype,~func[5],~func[4], func[3],~func[2],~func[1],~func[0]);
  and (i_srlv, rtype, ~func[5], ~func[4], ~func[3], func[2], func[1], ~func[0]);
  and (i_slt, rtype, func[5], ~func[4], func[3], ~func[2], func[1], ~func[0]);
  and (i_addi,~op[5],~op[4], op[3],~op[2],~op[1],~op[0]); // i format
  and (i_andi,~op[5],~op[4], op[3], op[2],~op[1],~op[0]);
  and (i_ori, ~op[5],~op[4], op[3], op[2],~op[1], op[0]);
  and (i_xori,~op[5],~op[4], op[3], op[2], op[1],~op[0]);
  and (i_lw, op[5],~op[4],~op[3],~op[2], op[1], op[0]);
  and (i_sw, op[5],~op[4], op[3],~op[2], op[1], op[0]);
  and (i_beq, ~op[5],~op[4],~op[3], op[2],~op[1],~op[0]);
  and (i_bne, ~op[5],~op[4],~op[3], op[2],~op[1], op[0]);
  and (i_lui, ~op[5],~op[4], op[3], op[2], op[1], op[0]);
  and (i_bgtz, ~op[5], ~op[4], ~op[3], op[2], op[1], op[0]);
  and (i_j, ~op[5],~op[4],~op[3],~op[2], op[1],~op[0]); // j format
  and (i_jal, ~op[5],~op[4],~op[3],~op[2], op[1], op[0]);
  
  logic i_rs, i_rt;
  // instructions that use rs
  assign i_rs = i_add | i_sub | i_and | i_or | i_xor | i_jr | i_addi |
              i_andi | i_ori | i_xori | i_lw | i_sw | i_beq | i_bne | i_bgtz;
  // instructions that use rt
  assign i_rt = i_add | i_sub | i_and | i_or | i_xor | i_sll | i_srl |
              i_sra | i_sw | i_beq | i_bne |i_srlv;
  // pipeline stall caused by data dependency with lw instruction
  assign nostall = ~(ewreg & em2reg & (ern != 0) & (i_rs & (ern == rs) | i_rt & (ern == rt)));

  logic [1:0] fwda, fwdb; // forwarding, multiplexer's select signals
  always_comb begin
    // forward control signal for alu input a
    fwda = 2'b00; // default: no hazards
    if (ewreg & (ern != 0) & (ern == rs) & ~em2reg) begin
      fwda = 2'b01; // select exe_alu
    end else begin
      if (mwreg & (mrn != 0) & (mrn == rs) & ~mm2reg) begin
        fwda = 2'b10; // select mem_alu
      end else begin
        if (mwreg & (mrn != 0) & (mrn == rs) & mm2reg) begin
          fwda = 2'b11; // select mem_lw
        end
      end
    end
    // forward control signal for alu input b
    fwdb = 2'b00; // default: no hazards  
    if (ewreg & (ern != 0) & (ern == rt) & ~em2reg) begin
      fwdb = 2'b01; // select exe_alu
    end else begin
      if (mwreg & (mrn != 0) & (mrn == rt) & ~mm2reg) begin
        fwdb = 2'b10; // select mem_alu
      end else begin
        if (mwreg & (mrn != 0) & (mrn == rt) & mm2reg) begin
          fwdb = 2'b11; // select mem_lw
        end
      end
    end
  end
  // control signals
  assign wreg =(i_add |i_sub |i_and |i_or |i_xor |i_sll |i_srl |
                i_sra |i_addi|i_andi|i_ori |i_xori|i_lw |i_lui |
                i_jal | i_srlv)& nostall; // prevent from executing twice
  assign regrt = i_addi|i_andi|i_ori |i_xori|i_lw |i_lui;
  assign jal = i_jal;
  assign m2reg = i_lw;
  assign shift = i_sll |i_srl |i_sra;
  assign aluimm = i_addi|i_andi|i_ori |i_xori|i_lw |i_lui |i_sw;
  assign sext = i_addi|i_lw |i_sw |i_beq |i_bne;
  assign aluc[4] = i_sra;
  assign aluc[3] = i_sub |i_or |i_srl |i_sra |i_ori |i_lui |i_srlv |i_slt;
  assign aluc[2] = i_slt;
  assign aluc[1] = i_xor |i_sll |i_srl |i_sra |i_xori|i_beq |i_bne |i_lui |i_srlv;
  assign aluc[0] = i_and |i_or |i_sll |i_srl |i_sra |i_andi |i_ori |i_srlv;
  assign wmem = i_sw & nostall; // prevent from executing twice
  assign pcsrc[1] = i_jr |i_j |i_jal;
  assign pcsrc[0] = i_beq & rsrtequ | i_bne & ~rsrtequ | i_j | i_jal | i_bgtz & bgtz_flag;
endmodule

module pipedereg (input logic dwreg, dm2reg, dwmem, input logic [4:0] daluc, input logic daluimm, input logic [31:0] da, db, dimm, input logic [4:0] drn, input logic dshift,
                  djal, input logic [31:0] dpc4, input logic clk, clrn, output logic ewreg, em2reg, ewmem, output logic [4:0] ealuc, output logic ealuimm, output logic [31:0] ea,
                  eb, eimm, output logic [4:0] ern, output logic eshift, ejal, output logic [31:0] epc4); // ID/EXE pipeline register

  always_ff @(negedge clrn or posedge clk) begin
    if (!clrn) begin // clear
      ewreg <= 0;   em2reg <= 0;
      ewmem <= 0;   ealuc <= 0;
      ealuimm <= 0; ea <= 0;
      eb <= 0;      eimm <= 0;
      ern <= 0;     eshift <= 0;
      ejal <= 0;    epc4 <= 0;
    end else begin // register
      ewreg <= dwreg;     em2reg <= dm2reg;
      ewmem <= dwmem;     ealuc <= daluc;
      ealuimm <= daluimm; ea <= da;
      eb <= db;           eimm <= dimm;
      ern <= drn;         eshift <= dshift;
      ejal <= djal;       epc4 <= dpc4;
    end
  end
endmodule

module pipeexe (input logic [4:0] ealuc, input logic ealuimm, input logic [31:0] ea, eb, eimm, input logic eshift, input logic [4:0] ern0, input logic [31:0] epc4, input logic ejal, output logic [4:0] ern, 
                output logic [31:0] ealu);

  logic [31:0] alua, alub, ealu0, epc8, esa; // alu input a // alu input b // alu result // pc+8 // shift amount
  logic z; // alu z flag, not used
  assign esa = {eimm[5:0], eimm[31:6]};
  cla32 ret_addr (epc8,epc4,32'h4,1'b0); // pc+8
  mux2 #(32) alu_in_a (ea,esa,eshift,alua); // alu input a
  mux2 #(32) alu_in_b (eb,eimm,ealuimm,alub); // alu input b
  mux2 #(32) save_pc8 (ealu0,epc8,ejal,ealu); // alu result or pc+8
  assign ern = ern0 | {5{ejal}}; // dest reg #, jal: 31
  alu al_unit (alua,alub,ealuc,ealu0,z); // alu result, z flag
endmodule

module pipeemreg (input logic ewreg, em2reg, ewmem, input logic [31:0] ealu, eb, input logic [4:0] ern, input logic clk, clrn, output logic mwreg, mm2reg,
                  mwmem, output logic [31:0] malu, mb, output logic [4:0] mrn); // EXE/MEM pipeline register

  always_ff @(negedge clrn or posedge clk) begin
    if (!clrn) begin // clear
      mwreg <= 0; mm2reg <= 0;
      mwmem <= 0; malu <= 0;
      mb <= 0; mrn <= 0;
    end else begin // register
      mwreg <= ewreg; mm2reg <= em2reg;
      mwmem <= ewmem; malu <= ealu;
      mb <= eb; mrn <= ern;
    end
  end
endmodule

module pipemem (input logic we, input logic [31:0] addr, datain, input logic clk, output logic [31:0] dataout); // MEM stage
  pl_data_mem dmem (clk,dataout,datain,addr,we); // data memory
endmodule

module pipemwreg (input logic mwreg, mm2reg, input logic [31:0] mmo, malu, input logic [4:0] mrn, input logic clk, clrn, output logic wwreg, wm2reg, output logic [31:0] wmo, walu,
                  output logic [4:0] wrn); // MEM/WB pipeline register

  always_ff @(negedge clrn or posedge clk) begin
    if (!clrn) begin // clear
      wwreg <= 0; wm2reg <= 0;
      wmo <= 0; walu <= 0;
      wrn <= 0;
    end else begin // register
      wwreg <= mwreg; wm2reg <= mm2reg;
      wmo <= mmo; walu <= malu;
      wrn <= mrn;
    end
  end
endmodule

module pipewb (input logic [31:0] walu, wmo, input logic wm2reg, output logic [31:0] wdi); // WB stage // alu result or pc+8 in WB stage // data out (from mem) in WB stage // memory to register in WB stage // data to be written into regfile
  mux2 #(32) wb (walu,wmo,wm2reg,wdi); // select for wdi
endmodule

module pl_inst_mem (input logic [31:0] a, output logic [31:0] inst); // instruction memory, rom

  logic [31:0] rom [0:63]; // rom cells: 64 words * 32 bits
  // rom[word_addr] = instruction // (pc) label instruction
  assign rom[6'h00] = 32'h20090008; // (00) main: lui $1, 0
  assign rom[6'h01] = 32'h200affff; // (04) ori $4, $1, 80
  assign rom[6'h02] = 32'h1D200002; // (08) call: jal sum
  assign rom[6'h03] = 32'h00000000; // (44) dslot2: nop
  assign rom[6'h04] = 32'h214a0014; // (0c) dslot1: addi $5, $0, 4
  assign rom[6'h05] = 32'h214afffb; // (10) return: sw $2, 0($4)
  assign rom[6'h06] = 32'h20060019; // (14) lw $9, 0($4)
  assign rom[6'h07] = 32'h20080001; // (18) sub $8, $9, $4
  assign rom[6'h08] = 32'h20020005; // (1c) addi $5, $0, 3
  assign rom[6'h09] = 32'h2003000c; // (20) loop2: addi $5, $5, -1
  assign rom[6'h0a] = 32'h2067fff7; // (24) ori $8, $5, 0xffff
  assign rom[6'h0b] = 32'h00e22025; // (28) xori $8, $8, 0x5555
  assign rom[6'h0c] = 32'h00642824; // (2c) addi $9, $0, -1
  assign rom[6'h0d] = 32'h00a42820; // (30) andi $10,$9,0xffff
  assign rom[6'h0e] = 32'h1d400013; // (34) or $6, $10, $9
  assign rom[6'h0f] = 32'h00000000; // (44) dslot2: nop
  assign rom[6'h10] = 32'h0064202a; // (38) xor $8, $10, $9
  assign rom[6'h11] = 32'h10800002; // (3c) and $7, $10, $6
  assign rom[6'h12] = 32'h00000000; // (44) dslot2: nop
  assign rom[6'h13] = 32'h20050000; // (40) beq $5, $0, shift
  assign rom[6'h14] = 32'h00e2202a; // (44) dslot2: nop
  assign rom[6'h15] = 32'h00853820; // (48) j loop2
  assign rom[6'h16] = 32'h00e23822; // (4c) dslot3: nop
  assign rom[6'h17] = 32'h00e63826; // (50) shift: addi $5, $0, -1
  assign rom[6'h18] = 32'h34e7000f; // (54) sll $8, $5, 15
  assign rom[6'h19] = 32'h30e7000f; // (58) sll $8, $8, 16
  assign rom[6'h1a] = 32'h38e7001e; // (5c) sra $8, $8, 16
  assign rom[6'h1b] = 32'h00073882; // (60) srl $8, $8, 15
  assign rom[6'h1c] = 32'h01073806; // (64) finish: j finish
  assign rom[6'h1d] = 32'hac670044; // (68) dslot4: nop
  assign rom[6'h1e] = 32'h8c020050; // (6c) sum: add $8, $0, $0
  assign rom[6'h1f] = 32'h20060084; // (70) loop: lw $9, 0($4)
  assign rom[6'h20] = 32'h00c00008; // (74) stall: add $8, $8, $9
  assign rom[6'h21] = 32'h20020001; // (78) addi $5, $5, -1
  assign rom[6'h22] = 32'hac020054; // (7c) bne $5, $0, loop
  assign rom[6'h23] = 32'h00000000; // (44) dslot2: nop
  //assign rom[6'h20] = 32'h20840004; // (80) dslot5: addi $4, $4, 4
  //assign rom[6'h21] = 32'h03e00008; // (84) jr $31
  //assign rom[6'h22] = 32'h00081000; // (88) dslot6: sll $2, $8, 0
  assign inst = rom[a[7:2]]; // use 6-bit word address to read rom
endmodule

module pl_data_mem (input logic clk, output logic [31:0] dataout, input logic [31:0] datain, addr, input logic we); // data memory, ram

  logic [31:0] ram [0:63]; // ram cells: 32 words * 32 bits
  assign dataout = ram[addr[31:2]]; // use 5-bit word address
  always_ff @(posedge clk) begin
    if (we) ram[addr[31:2]] = datain; // write ram
  end
  integer i;
  initial begin // ram initialization
    for (i = 0; i < 64; i = i + 1)
      ram[i] = 0;
    // ram[word_addr] = data // (byte_addr) item in data array
    //ram[5'h14] = 32'h000000a3; // (50) data[0] 0 + a3 = a3
    //ram[5'h15] = 32'h00000027; // (54) data[1] a3 + 27 = ca
    //ram[5'h16] = 32'h00000079; // (58) data[2] ca + 79 = 143
    //ram[5'h17] = 32'h00000115; // (5c) data[3] 143 + 115 = 258
    // ram[5'h18] should be 0x00000258, the sum stored by sw instruction
  end
endmodule


module alu (input logic [31:0] a, b, input logic [4:0] aluc, output logic [31:0] r, output logic z); // 32-bit alu with a zero flag

  logic [31:0] d_and, d_or, d_xor, d_lui, d_and_or, d_xor_lui, d_as, d_sh;
  assign d_and = a & b; // x 0 0 1 AND
  assign d_or = a | b; // x 1 0 1 OR
  assign d_xor = a ^ b; // x 0 1 0 XOR
  assign d_lui = {b[15:0],16'h0}; // x 1 1 0 LUI
  assign d_and_or = aluc[3]? d_or : d_and; // 0 0 1 1 SLL
  assign d_xor_lui = aluc[3]? d_lui : d_xor; // 0 1 1 1 SRL

  // addsub32 (a,b,sub, s);
  addsub32 as32 (d_as,a,b,aluc[3]); // add/sub
  // shift (d,sa, right, arith, sh);
  shift shifter (b,a[4:0],aluc[3],aluc[4],d_sh); // shift
  // mux4x32 (a0, a1, a2, a3, s, y);
  mux8 #(32) res (d_as, d_and_or, d_xor_lui, d_sh, d_as[31],32'dx, 32'dx, 32'dx, aluc[2:0], r); // alu result
  assign z = ~|r; // z = (r == 0)
endmodule

module shift (input logic [31:0] d, input logic [4:0] sa, input logic right, arith, output logic [31:0] sh); // barrel shift, behavioral style

  always_comb begin // always block
    if (!right) begin // if shift left
      sh = d << sa; // shift left sa bits
    end else if (!arith) begin // if shift right logical
      sh = d >> sa; // shift right logical sa bits
    end else begin // if shift right arithmetic
      sh = $signed(d) >>> sa; // shift right arithmetic sa bits
    end
  end
endmodule

module regfile (input logic [4:0] rna, rnb, input logic [31:0] d, input logic [4:0] wn, input logic we, clk, clrn, output logic [31:0] qa, qb); // 32x32 regfile

  logic [31:0] register [1:31]; // 31 32-bit registers
  assign qa = (rna == 0)? 0 : register[rna]; // read port A
  assign qb = (rnb == 0)? 0 : register[rnb]; // read port B
  integer i;
  always_ff @(posedge clk or negedge clrn) // write port
    if (!clrn) begin
      for (i = 1; i < 32; i = i + 1)
        register[i] <= 0; // reset
    end
    else begin
      if ((wn != 0) && we) // not reg[0] & enabled
        register[wn] <= d; // write d to reg[wn]
    end
endmodule

module mux2 #(parameter WIDTH = 8)
             (input  logic [WIDTH-1:0] d0, d1, 
              input  logic             s, 
              output logic [WIDTH-1:0] y);

  assign y = s ? d1 : d0; 
endmodule

module mux4 #(parameter WIDTH = 8)
             (input  logic [WIDTH-1:0] d0, d1, d2, d3, 
              input  logic [1:0]       s, 
              output logic [WIDTH-1:0] y);

  always_comb begin
    case(s)
        2'b00: y = d0;
        2'b01: y = d1;
        2'b10: y = d2;
        2'b11: y = d3;
    endcase
  end 
endmodule

module mux8 #(parameter WIDTH = 8)
             (input  logic [WIDTH-1:0] d0, d1, d2, d3, d4, d5, d6, d7, 
              input  logic [2:0]       s, 
              output logic [WIDTH-1:0] y);

  always_comb begin
    case(s)
        3'b000: y = d0;
        3'b001: y = d1;
        3'b010: y = d2;
        3'b011: y = d3;
        3'b100: y = d4;
        3'b101: y = d5;
        3'b110: y = d6;
        3'b111: y = d7;
    endcase
  end 
endmodule

module dffe32 (input logic [31:0] d, input logic clk, clrn, e, output logic [31:0] q); // a 32-bit register

  always_ff @(negedge clrn or posedge clk) begin
    if (!clrn) q <= 0; // q = 0 if reset
    else if (e) q <= d; // save d if enabled
  end
endmodule