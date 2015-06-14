canonical                    =  0.690454802741548
RK2_noadjust_final_error_list=[ 0.219892345371736 , 0.512249663253732 , 0.590530530954006 , 0.654821367228458 , 0.674945224049668 , 0.681633292809444 , 0.685788608102907 , 0.688252722091707]
RK2_adjusted_final_error_list=[ 0.665132847685988 , 0.686412554557701 , 0.689496546475445 , 0.690218056375439 , 0.690394681697134 , 0.690439677806910 , 0.690451078150846 , 0.690453861228918]
RK4_noadjust_final_error_list=[ 0.217663942041021 , 0.512835849305701 , 0.590653351650759 , 0.654844605046700 , 0.674950508898182 , 0.681634575948395 , 0.685788923240616 , 0.688252800061849]
RK4_adjusted_final_error_list=[ 0.685646521581585 , 0.690343683676600 , 0.690453100431482 , 0.690454780486604 , 0.690454802850925 , 0.690454803198193 , 0.690454803198701 , 0.690454803221983]


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
