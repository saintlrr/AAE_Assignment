# AAE_Assignment
Step 1.  We initialize const variables and read the rcvr.dat and eph.dat in the tables and add the table name lists in order. 
 
Step 2. Calculate the Position of each satellite

	Calculate the approximately satellite clock offset without consideration of relativistic effect according to ICD.pdf, Page 106, eq(2)
Δt_sv=a_(f_0 )+a_f1 (t-t_oc )+a_(f_2 ) (t-t_oc )^2
 
	Merge two tables
 
	Because the pseudorange is from the satellite at transmit the signal time to receiver at receive time, we should calculate the satellites’ positions at transmission time.
 
	To specific to calculate the satellites’ positions based on ICD.pdf page 116-118
 
We can get the function to calculate satellite positions
 
Finally, using above function we can get satellites’ positions[x,y,z] at transmission time based on the svid order[5,6,17,30,10,23,22,26]
-8855451.25241948	-22060172.8092243	-11922096.6034882
-8087137.79683860	-16946009.4227108	18816191.7549088
-21277076.7306835	-7467236.05896071	14287508.1869877
-17713787.7954120	-19797565.6371920	19203.6500188289
9027722.46659624	-12319176.0234634	21737388.2917114
-19452220.1443220	-16750493.6910820	-6918515.43486897
-13649579.4737585	8229423.98791781	21122957.9604718
6163064.79386048	-25286737.6137436	-3541185.66549232

	Recalculate the satellite clock offset with relativistic correction based on ICD.pdf page 106
 
We can get the satellites’ clock offset based on the svid order[5,6,17,30,10,23,22,26]
0.000189065737037372
-8.39324641464491e-08
-0.000204902643239992
-1.00411322004581e-05
3.32476279119490e-05
1.03602104492918e-05
0.000222678092585821
0.000280993171239945

	Calculate the initial iteration δx
δρ=ρ+c⋅δt^sat-ρ ̂-c⋅δt_rcvr
H=[(x_i-χ^sat)/ρ ̂ ,(y_i-y^sat)/ρ ̂ ,(z_i-z^sat)/ρ ̂ ,1]
δχ=(H^T H)^(-1) H^T δρ
 
Thus, we can get the δx_0
-5714.33483899312
1080.96401795638
-2605.92555469761
519462.246977505

	Do the iterations until δx<10^(-4)
 
Finally, we can get the receiver position and receiver clock offset
-2700399.93620375
-4292561.50063647
3855273.10354145
519461.693994075
