module mccomp (input logic clk,clrn, output logic [2:0] q, output logic [31:0] a, b, alu, adr, tom, fromm, pc, ir, input logic memclk); // clock and reset memory  // state  // alu input a // alu input b // alu result //  memory address // data to memory// data from memory // program counter // instruction register // clk for sync ram;
    logic wmem; // memory write enable
    mccpu mc_cpu (clk,clrn,fromm,pc,ir,a,b,alu,wmem,adr,tom,q); // cpu
    iordmem memory (memclk,wmem,adr,tom,fromm); // memory
endmodule

module mccpu (input logic clk,clrn, input logic [31:0] frommem, output logic [31:0] pc,inst,alua,alub,alur,output logic wmem, output logic [31:0] madr,tomem, output logic [2:0] state);
              // clock and reset // data from memory // program counter // instruction // alu input a // alu input b // alu result // memory write enable // memory address // data to memory // state

// instruction fields
    logic [5:0] op;
    logic [4:0] rs;
    logic [4:0] rt;
    logic [4:0] rd;
    logic [5:0] func;
    logic [15:0] imm;
    logic [25:0] addr;

    assign op = inst[31:26]; // op
    assign rs = inst[25:21]; // rs
    assign rt = inst[20:16]; // rt
    assign rd = inst[15:11]; // rd
    assign func = inst[05:00]; // func
    assign imm = inst[15:00]; // immediate
    assign addr = inst[25:00]; // address
// control signals
    logic [4:0] aluc; // alu operation control
    logic [1:0] pcsrc; // select pc source
    logic wreg; // write regfile
    logic regrt; // dest reg number is rt
    logic m2reg; // instruction is an lw
    logic shift; // instruction is a shift
    logic [1:0] alusrcb; // alu input b selection
    logic jal; // instruction is a jal
    logic sext; // is sign extension
    logic wpc; // write pc
    logic wir; // write ir
    logic iord; // select memory address
    logic selpc; // select pc
// datapath wires
    logic [31:0] bpc; // branch target address
    logic [31:0] npc; // next pc
    logic [31:0] qa; // regfile output port a
    logic [31:0] qb; // regfile output port b
    logic [31:0] wd; // regfile write port data
    logic [31:0] r; // alu out or mem
    logic [31:0] sa;
    logic [15:0] s16;
    logic [31:0] i32;
    logic [31:0] dis;
    logic [31:0] jpc;
    logic [4:0] reg_dest; // rs or rt
    logic [4:0] wn;
    logic z; // alu, zero tag
    logic bgtz_flag; //bgtz, bgtz flag
    logic [31:0] rega; // register a
    logic [31:0] regb; // register b
    logic [31:0] regc; // register c
    logic [31:0] data; // output of dr
    logic [31:0] opa; // sa or output of reg a

    assign sa = {27'b0,inst[10:6]}; // shift amount
    assign s16 = {16{sext & inst[15]}}; // 16-bit signs
    assign i32 = {s16,imm}; // 32-bit immediate
    assign dis = {s16[13:0],imm,2'b00}; // word distance
    assign jpc = {pc[31:28],addr,2'b00}; // jump target address
    assign wn = reg_dest | {5{jal}}; // regfile write reg #
    // datapath
    mccu control_unit (op,func,z,bgtz_flag,clk,clrn, // control unit
                      wpc,wir,wmem,wreg,
                      iord,regrt,m2reg,aluc,
                      shift,selpc,alusrcb,
                      pcsrc,jal,sext,state);
    dffe32 ip (npc,clk,clrn,wpc,pc); // pc register
    dffe32 ir (frommem,clk,clrn,wir,inst); // instruction register
    dff32 dr (frommem,clk,clrn,data); // data register
    dff32 reg_a (qa, clk,clrn,rega); // register a
    dff32 reg_b (qb, clk,clrn,regb); // register b
    dff32 reg_c (alur,clk,clrn,regc); // register c
    mux2 #(32) aorsa (rega,sa,shift,opa); // sa or output of reg a
    mux2 #(32) alu_a (opa,pc,selpc,alua); // alu input a
    mux4 #(32) alu_b (regb,32'h4,i32,dis,alusrcb,alub); // alu input b
    mux2 #(32) alu_m (regc,data,m2reg,r); // alu out or mem data
    mux2 #(32) mem_a (pc,regc,iord,madr); // memory address
    mux2 #(32) link (r,pc,jal,wd); // r or pc
    mux2 #(5) reg_wn (rd,rt,regrt,reg_dest); // rs or rt
    mux4 #(32) nextpc(alur,regc,qa,jpc,pcsrc,npc); // next pc
    regfile rf (rs,rt,wd,wn,wreg,clk,clrn,qa,qb); // register file
    alu alunit (alua,alub,aluc,alur,z,bgtz_flag); // alu
    assign tomem = regb; // output of reg b
endmodule

module mccu (input logic [5:0] op,func, input logic z,bgtz_flag,clk,clrn, output logic wpc,wir,wmem,wreg,iord,regrt,m2reg, output logic [4:0] aluc,
            output logic shift,alusrca, output logic [1:0] alusrcb,pcsrc, output logic jal,sext, output logic [2:0] state); // control unit
// op, func // alu zero flag // clock and reset // write pc // write ir // memory write enable // write regfile // select memory address // dest reg number is rt // instruction is an lw // alu operation control
    logic [2:0] next_state; // next state
    parameter [2:0] sif = 3'b000, // IF state
    sid = 3'b001, // ID state
    sexe = 3'b010, // EXE state
    smem = 3'b011, // MEM state
    swb = 3'b100; // WB state
    logic rtype,i_add,i_sub,i_and,i_or,i_xor,i_sll,i_srl,i_sra,i_jr, i_srlv, i_slt;
    logic i_addi,i_andi,i_ori,i_xori,i_lw,i_sw,i_beq,i_bne,i_lui,i_j,i_jal, i_bgtz;
    // and (out, in1, in2, ...); // instruction decode
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
    and (i_srlv, rtype,~func[5],~func[4], ~func[3],func[2],func[1],~func[0]);
    and (i_slt, rtype,func[5],~func[4], func[3],~func[2],func[1],~func[0]);
    and (i_addi,~op[5],~op[4], op[3],~op[2],~op[1],~op[0]); // i format
    and (i_andi,~op[5],~op[4], op[3], op[2],~op[1],~op[0]);
    and (i_ori, ~op[5],~op[4], op[3], op[2],~op[1], op[0]);
    and (i_xori,~op[5],~op[4], op[3], op[2], op[1],~op[0]);
    and (i_lw, op[5],~op[4],~op[3],~op[2], op[1], op[0]);
    and (i_sw, op[5],~op[4], op[3],~op[2], op[1], op[0]);
    and (i_beq, ~op[5],~op[4],~op[3], op[2],~op[1],~op[0]);
    and (i_bne, ~op[5],~op[4],~op[3], op[2],~op[1], op[0]);
    and (i_lui, ~op[5],~op[4], op[3], op[2], op[1], op[0]);
    and (i_bgtz,~op[5],~op[4], ~op[3], op[2], op[1], op[0]);
    and (i_j, ~op[5],~op[4],~op[3],~op[2], op[1],~op[0]); // j format
    and (i_jal, ~op[5],~op[4],~op[3],~op[2], op[1], op[0]);
    logic i_shift = i_sll | i_srl | i_sra;
    always_comb begin // default outputs:
        wpc = 0; // do not write pc
        wir = 0; // do not write ir
        wmem = 0; // do not write memory
        wreg = 0; // do not write regfile
        iord = 0; // select pc as address
        aluc = 5'bx0000; // alu operation: add
        alusrca = 0; // alu a: reg a or sa
        alusrcb = 2'h0; // alu input b: reg b
        regrt = 0; // reg dest no: rd
        m2reg = 0; // select reg c
        shift = 0; // select reg a
        pcsrc = 2'h0; // select alu output
        jal = 0; // not a jal
        sext = 1; // sign extend
        case (state) //------------------------------------------------- IF:
            sif: begin // IF state
                wpc = 1; // write PC
                wir = 1; // write IR
                alusrca = 1; // PC
                alusrcb = 2'h1; // 4
                next_state = sid; // next state: ID
            end //------------------------------------------------------ ID:
            sid: begin // ID state
                if (i_j) begin // j instruction
                    pcsrc = 2'h3; // jump address
                    wpc = 1; // write PC
                    next_state = sif; // next state: IF
                end else if (i_jal) begin // jal instruction
                    pcsrc = 2'h3; // jump address
                    wpc = 1; // write PC
                    jal = 1; // reg no = 31
                    wreg = 1; // save PC+4
                    next_state = sif; // next state: IF
                end else if (i_jr) begin // jr instruction
                    pcsrc = 2'h2; // jump register
                    wpc = 1; // write PC
                    next_state = sif; // next state: IF
                end else begin // other instructions
                    aluc = 5'bx0000; // add
                    alusrca = 1; // PC
                    alusrcb = 2'h3; // branch offset
                    next_state = sexe; // next state: EXE
                end
            end //----------------------------------------------------- EXE:
            sexe: begin // EXE state
                aluc[4] = i_sra;
                aluc[3] = i_sub |i_or |i_srl |i_sra |i_ori |i_lui |i_srlv |i_slt;
                aluc[2] = i_slt;
                aluc[1] = i_xor |i_sll |i_srl |i_sra |i_xori|i_beq |i_bne |i_lui |i_srlv;
                aluc[0] = i_and |i_or |i_sll |i_srl |i_sra |i_andi |i_ori |i_srlv;
                if (i_beq || i_bne || i_bgtz) begin // beq or bne inst
                    pcsrc = 2'h1; // branch address
                    wpc = i_beq & z | i_bne & ~z | i_bgtz & bgtz_flag; // write PC
                    next_state = sif; // next state: IF
                end else begin // other instruction
                    if (i_lw || i_sw) begin // lw or sw inst
                        alusrcb = 2'h2; // select offset
                        next_state = smem; // next state: MEM
                    end else begin // other instruction
                    if (i_shift) shift = 1; // shift instruction
                    if (i_addi || i_andi || i_ori || i_xori || i_lui) alusrcb = 2'h2; // select immediate
                    if (i_andi || i_ori || i_xori) sext=0; // 0 extend
                        next_state = swb; // next state: WB
                    end
                end
            end //----------------------------------------------------- MEM:
            smem: begin // MEM state
                iord = 1; // memory address = C
                if (i_lw) begin // load
                    next_state = swb; // next state: WB
                end else begin // store
                    wmem = 1; // write memory
                    next_state = sif; // next state: IF
                end
            end //------------------------------------------------------ WB:
            swb: begin // WB state
                if (i_lw) m2reg = 1; // select memory data
                if (i_lw || i_addi || i_andi || i_ori || i_xori || i_lui)
                    regrt = 1; // reg dest no: rt
                wreg = 1; // write register file
                next_state = sif; // next state: IF
            end //------------------------------------------------------ END
            default: begin
                next_state = sif; // default state: IF
            end
        endcase
    end
    always_ff @(posedge clk or negedge clrn) begin
        if (!clrn) begin
            state <= sif; // reset state to IF
        end else begin
            state <= next_state; // state transition
        end
    end
endmodule

module alu (input logic [31:0] a, b, input logic [4:0] aluc, output logic [31:0] r, output logic z, bgtz_flag); // 32-bit alu with a zero flag

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
  mux8 #(32) res (d_as, d_and_or, d_xor_lui, d_sh, d_as[31],32'dz, 32'dz, 32'dz, aluc[2:0], r); // alu result
  assign z = ~|r; // z = (r == 0)
  assign bgtz_flag = (!(a[31] == 1'b1)) && (a > 32'b0); //to check if greater than zero and also not zero
endmodule

/*
module alu (input logic [31:0] a, b, input logic [3:0] aluc, output logic [31:0] r, output logic z); // 32-bit alu with a zero flag

  logic [31:0] d_and, d_or, d_xor, d_lui, d_and_or, d_xor_lui, d_as, d_sh;
  assign d_and = a & b; // x 0 0 1 AND
  assign d_or = a | b; // x 1 0 1 OR
  assign d_xor = a ^ b; // x 0 1 0 XOR
  assign d_lui = {b[15:0],16'h0}; // x 1 1 0 LUI
  assign d_and_or = aluc[2]? d_or : d_and; // 0 0 1 1 SLL
  assign d_xor_lui = aluc[2]? d_lui : d_xor; // 0 1 1 1 SRL

  // addsub32 (a,b,sub, s);
  addsub32 as32 (d_as,a,b,aluc[2]); // add/sub
  // shift (d,sa, right, arith, sh);
  shift shifter (b,a[4:0],aluc[2],aluc[3],d_sh); // shift
  // mux4x32 (a0, a1, a2, a3, s, y);
  mux4 #(32) res (d_as, d_and_or, d_xor_lui, d_sh, aluc[1:0], r); // alu result
  assign z = ~|r; // z = (r == 0)
endmodule
*/

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

module dff32 (input logic [31:0] d, input logic clk, clrn, output logic [31:0] q); // a 32-bit register

  always_ff @(negedge clrn or posedge clk) begin
    if (!clrn) q <= 0; // q = 0 if reset
    else q <= d; // save d if enabled
  end
endmodule

module iordmem (input  logic    clk, we,
            input  logic [31:0]  a, 
            input  logic [31:0] wd,
            output logic [31:0] rd);
  
  logic [31:0] RAM[0:63];

  initial begin
      $readmemh("D:/Downloads/MIPS/memfile.mem",RAM);
  end
  
  always_ff @(posedge clk) begin
    if(we) RAM[a[7:2]] <= wd;
  end

  assign rd = RAM[a[7:2]];
endmodule