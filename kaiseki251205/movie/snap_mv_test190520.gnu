#set terminal postscript enhanced color eps "Times-ltalic" 20 
set terminal png enhanced
#m0からm1までのファイルをm2毎に連番処理する
m0=1181
m1=1200
m2=20

###BOX###
set size ratio -1
set pm3d explicit map
set xlabel "R(R_S)"
set ylabel "z(R_S)"
set xrange [0:750000]
set yrange [0:750000]
#set xtics(-1000,-500,0,500,1000)
#set ytics(-1000,-500,0,500,1000)
unset key 
#########

#ベクトルのスタイル
set style arrow 1 nofilled linewidth 1.5 lc rgb 'black' #データ
set style arrow 2 nofilled linewidth 1.5 lc rgb 'black' #見本

#旧式if (exist("n")==0 || n<0) n=0 ;
#ループでファイルを順番に処理する
do for [n=m0:m1:m2] {

   file0=sprintf("mvdata_%.04d.dat",n)
   file1=sprintf("vecdata_%.04d.dat",n)
   #outfile1=sprintf("density_%.04d.eps",n/20.)
   outfile1=sprintf("rho_No.2_v_%.04d.png",n) #(/40.0)
   #outfile1=sprintf("rho_22-18_%.04d.png",n) #(/40.0)
   outfile2=sprintf("temp_No.2_%.04d.png",n) #(/40.0)
   outfile3=sprintf("vr_No.2_%.04d.png",n) 
   outfile4=sprintf("Fradr_No.2_%.04d.png",n) 
   outfile5=sprintf("forceM_No.2_%.04d.png",n) 
   


   set format cb "%T"  
   set log z
   set log cb
   
   #温度プロット
   set palette defined ( 0 '#352a87',1 '#2053d4',2 '#0d75dc',3 '#0c93d2',4 '#07a9c2',5 '#38b99e',6 '#7cbf7b',7 '#b7bd64',8 '#f1b94a',9 '#fad32a',10 '#f9fb0e')
   set cbrange[1.e4:1.e+9]
   #set cbrange[-0.5:1.0]
   
   set output outfile2
   splot file0 us 1:2:4 with pm3d


   #密度プロット
   set palette define (1 "blue", 2 "cyan", 3 "yellow", 4 "orange", 5 "red")
   #set zrange [1e-23:1e-9]
   set cbrange [1.e-22:1.e-18]

   set output outfile1
   #set label 2 "log {/Symbol r}(g cm^{-3})" at 2250,0 center
   #set label 2 rotate by 90
   set arrow 1 from 7*1.e5,6*1.e5 to 7*1.e5,6.75*1.e5 arrowstyle 2
   set label 1 "0.01 c" at 5.8*1.e5,6.3*1.e5
   splot file0 us 1:2:3 with pm3d,file1 us 1:2:($3+1.e-17):($4*7.5*1.e6):($5*7.5*1.e6):($6*7.5*1.e6) with vector arrowstyle 1
   #,\


   #vrプロット
   set palette define (1 "blue", 2 "white", 3 "red")
   #set zrange [1e-23:1e-5]
   
   set format cb "%g" 
   unset log z
   unset log cb
   set cbrange [-0.02:0.02]
   set output outfile3
   splot file0 us 1:2:5 with pm3d

   #Frad_r/grav_rプロット
   set log cb
   set palette define (1 "blue", 2 "cyan", 3 "yellow", 4 "orange", 5 "red")
   #set zrange [1e-23:1e-9]
   set cbrange [1.e-1:100]

   set output outfile4
   splot file0 us 1:2:6 with pm3d


   #force_Mプロット
   set palette define (1 "blue", 2 "cyan", 3 "yellow", 4 "orange", 5 "red")
   #set zrange [1e-23:1e-9]
   #set cbrange [0.0:]
   set cbrange [1.e-5:0.1]

   set output outfile5
   splot file0 us 1:2:7 with pm3d
   unset log cb

   

}

#if (n<m1) n=n+m2; reread
