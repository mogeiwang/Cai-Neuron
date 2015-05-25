
canonical                     =0.175578230757126
RK2_noadjust_final_error_list=[ 0.905871657082677 , 0.982679406940519 , 0.048688959563439 , 0.108893854748690 , 0.145775348112375 , 0.161937648514800 , 0.169672260781170 , 0.172953565766778]
RK2_adjusted_final_error_list=[ 0.156209665208635 , 0.170908792882701 , 0.174424795225142 , 0.175291832140348 , 0.175506118784725 , 0.175560285433530 , 0.175573767461711 , 0.175577118802827]
RK4_noadjust_final_error_list=[ 0.906503285509976 , 0.982822531786080 , 0.071995251519531 , 0.108903143689566 , 0.145777722981026 , 0.161938246801384 , 0.169672410923413 , 0.172953603358968]
RK4_adjusted_final_error_list=[ 0.175573571928253 , 0.175577947786528 , 0.175578214181014 , 0.175578230475127 , 0.175578231561049 , 0.175578231612355 , 0.175578231688056 , 0.175578231673770]



for i in range(len(RK2_noadjust_final_error_list)):
    RK2_noadjust_final_error_list[i] = abs(RK2_noadjust_final_error_list[i] - canonical)
    if (RK2_noadjust_final_error_list[i]>0.5): RK2_noadjust_final_error_list[i] = 1 - RK2_noadjust_final_error_list[i]

for i in range(len(RK2_adjusted_final_error_list)):
    RK2_adjusted_final_error_list[i] = abs(RK2_adjusted_final_error_list[i] - canonical)
    if (RK2_adjusted_final_error_list[i]>0.5): RK2_adjusted_final_error_list[i] = 1 - RK2_adjusted_final_error_list[i]

for i in range(len(RK4_noadjust_final_error_list)):
    RK4_noadjust_final_error_list[i] = abs(RK4_noadjust_final_error_list[i] - canonical)
    if (RK4_noadjust_final_error_list[i]>0.5): RK4_noadjust_final_error_list[i] = 1 - RK4_noadjust_final_error_list[i]

for i in range(len(RK4_adjusted_final_error_list)):
    RK4_adjusted_final_error_list[i] = abs(RK4_adjusted_final_error_list[i] - canonical)
    if (RK4_adjusted_final_error_list[i]>0.5): RK4_adjusted_final_error_list[i] = 1 - RK4_adjusted_final_error_list[i]


ts_list=[]
for i in range(8): ts_list.append(0.001 / 2**i)

figure()
plot(ts_list, RK2_noadjust_final_error_list, '+b', label='RK2 Org')
plot(ts_list, RK2_adjusted_final_error_list, '.b', label='RK2 Adj')
plot(ts_list, RK4_noadjust_final_error_list, 'xr', label='RK4 Org')
plot(ts_list, RK4_adjusted_final_error_list, '*r', label='RK4 Adj')

# ----------------------------------------------------------------
# time step lengthes: [0.001, 0.0005, 0.00025, 0.000125, 6.25e-05, 3.125e-05, 1.5625e-05]

for i in range(8):
    ts_list[i]=log10(ts_list[i])
    RK2_noadjust_final_error_list[i]=log10(RK2_noadjust_final_error_list[i])
    RK2_adjusted_final_error_list[i]=log10(RK2_adjusted_final_error_list[i])
    RK4_noadjust_final_error_list[i]=log10(RK4_noadjust_final_error_list[i])
    RK4_adjusted_final_error_list[i]=log10(RK4_adjusted_final_error_list[i])

figure()
plot(ts_list, RK2_noadjust_final_error_list, '+b', label='RK2 Org')
plot(ts_list, RK2_adjusted_final_error_list, '.b', label='RK2 Adj')
plot(ts_list, RK4_noadjust_final_error_list, 'xr', label='RK4 Org')
plot(ts_list, RK4_adjusted_final_error_list, '*r', label='RK4 Adj')

z2n = polyfit(ts_list, RK2_noadjust_final_error_list, 1); p2n = poly1d(z2n) # print(z1) [ 0.31004412  6.65525   ];;; print(p1) 0.31 x + 6.655
z2a = polyfit(ts_list, RK2_adjusted_final_error_list, 1); p2a = poly1d(z2a)
z4n = polyfit(ts_list, RK4_noadjust_final_error_list, 1); p4n = poly1d(z4n)
z4a = polyfit(ts_list[:5], RK4_adjusted_final_error_list[:5], 1); p4a = poly1d(z4a)
#
print "fit RK2 noAdjust:", p2n
print "fit RK2 Adjusted:", p2a
print "fit RK4 noAdjust:", p4n
print "fit RK4 Adjusted:", p4a

plot(ts_list, p2n(ts_list), 'b', label='%s'%p2n)
plot(ts_list, p2a(ts_list), 'b', label='%s'%p2a)
plot(ts_list, p4n(ts_list), 'r', label='%s'%p4n)
plot(ts_list[:5], p4a(ts_list[:5]), 'r', label='%s'%p4a)

legend(loc='upper left')
xlabel("$\log_{10}$(time_step)")
ylabel("$\log_{10}$(error)")
# title("")
savefig('final_error_list.pdf', dpi=600)

