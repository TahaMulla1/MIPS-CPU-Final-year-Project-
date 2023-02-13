module add(output logic g, p, s, input logic a, b, c);
    assign s = a ^ b ^ c;
    assign g = a & b;
    assign p = a | b; 
endmodule

module gp(output logic g_out, p_out, c_out, input logic [1:0] g, p, input logic c_in);
    assign g_out = g[1] | p[1] & g[0];
    assign p_out = p[1] & p[0];
    assign c_out = g[0] | p[0] & c_in;
endmodule

module cla_2(output logic [1:0] s, 
        output logic g_out, p_out, 
        input logic [1:0] a, b, 
        input logic c_in);

    logic [1:0] g,p;
    logic c_out;
    add a0 (g[0], p[0], s[0], a[0], b[0], c_in);
    add a1 (g[1], p[1], s[1], a[1], b[1], c_out);
    
    gp gp0 (g_out, p_out, c_out, g, p, c_in);
            
endmodule

module cla_4(output logic [3:0] s, 
        output logic g_out, p_out, 
        input logic [3:0] a, b, 
        input logic c_in);

    logic [1:0] g,p;
    logic c_out;
    cla_2 a0 (s[1:0],g[0],p[0],a[1:0],b[1:0],c_in);
    cla_2 a1 (s[3:2],g[1],p[1],a[3:2],b[3:2],c_out);
    
    gp gp0 (g_out, p_out, c_out, g, p, c_in);
            
endmodule

module cla_8(output logic [7:0] s, 
        output logic g_out, p_out, 
        input logic [7:0] a, b, 
        input logic c_in);

    logic [1:0] g,p;
    logic c_out;
    cla_4 a0 (s[3:0],g[0],p[0],a[3:0],b[3:0],c_in);
    cla_4 a1 (s[7:4],g[1],p[1],a[7:4],b[7:4],c_out);
    
    gp gp0 (g_out, p_out, c_out, g, p, c_in);
            
endmodule

module cla_16(output logic [15:0] s, 
        output logic g_out, p_out, 
        input logic [15:0] a, b, 
        input logic c_in);

    logic [1:0] g,p;
    logic c_out;
    cla_8 a0 (s[7:0],g[0],p[0],a[7:0],b[7:0],c_in);
    cla_8 a1 (s[15:8],g[1],p[1],a[15:8],b[15:8],c_out);
    
    gp gp0 (g_out, p_out, c_out, g, p, c_in);
            
endmodule

module cla_32(output logic [31:0] s, 
        output logic g_out, p_out, 
        input logic [31:0] a, b, 
        input logic c_in);

    logic [1:0] g,p;
    logic c_out;
    cla_16 a0 (s[15:0],g[0],p[0],a[15:0],b[15:0],c_in);
    cla_16 a1 (s[31:16],g[1],p[1],a[31:16],b[31:16],c_out);
    
    gp gp0 (g_out, p_out, c_out, g, p, c_in);
            
endmodule

module cla32(output logic[31:0] s, input logic [31:0] a, b, input logic ci);
    logic g_out, p_out;
    cla_32 cla (s, g_out, p_out, a, b, ci);
endmodule

module addsub32 (output logic [31:0] s, input logic [31:0] a,b, input logic sub);
    cla32 as32 (s, a, b, sub);
endmodule