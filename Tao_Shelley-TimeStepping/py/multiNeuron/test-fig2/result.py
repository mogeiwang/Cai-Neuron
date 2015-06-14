# /r/T/C/T/p/m/test-fig2> cat output.txt | grep VOL0
canonical                    =  0.175578230757
RK2_noadjust_final_error_list=[ 0.905871657083 , 0.982679406941 , 0.0486889595634 , 0.108893854749 , 0.145775348112 , 0.161937648515 , 0.169672260781 , 0.172953565767]
RK2_adjusted_final_error_list=[ 0.156209665209 , 0.170908792883 , 0.174424795225 , 0.17529183214 , 0.175506118785 , 0.175560285434 , 0.175573767462 , 0.175577118803]
RK4_noadjust_final_error_list=[ 0.90650328551 , 0.982822531786 , 0.0719952515195 , 0.10890314369 , 0.145777722981 , 0.161938246801 , 0.169672410923 , 0.172953603359]
RK4_adjusted_final_error_list=[ 0.175573571928 , 0.175577947787 , 0.175578214181 , 0.175578230475 , 0.175578231561 , 0.175578231612 , 0.175578231688 , 0.175578231674]

# /r/T/C/T/p/m/test-fig2> cat output.txt | grep VOL1
#  canonical                    =  0.60875477235
#  RK2_noadjust_final_error_list=[ 0.0959259373898 , 0.777855330175 , 0.129331488674 , 0.421703588932 , 0.516507203813 , 0.559313875081 , 0.582440044365 , 0.597539031475]
#  RK2_adjusted_final_error_list=[ 0.540899483655 , 0.59283112151 , 0.604863444207 , 0.607790834164 , 0.608514731786 , 0.608694865521 , 0.608739812146 , 0.608751033018]
#  RK4_noadjust_final_error_list=[ 0.250044560082 , 0.813987012379 , 0.169767161226 , 0.421713171449 , 0.516509623835 , 0.559314480749 , 0.58244019559 , 0.597539069238]
#  RK4_adjusted_final_error_list=[ 0.608737079016 , 0.608753699793 , 0.608754706077 , 0.608754767896 , 0.608754771812 , 0.608754771968 , 0.608754772093 , 0.608754771998]

# /r/T/C/T/p/m/test-fig2> cat output.txt | grep VOL2
#  canonical                    =  -0.958612424639
#  RK2_noadjust_final_error_list=[ -0.958613029344 , -0.958612574286 , -0.958612461896 , -0.958612433968 , -0.958612427006 , -0.958612425271 , -0.958612424832 , -0.958612424732]
#  RK2_adjusted_final_error_list=[ -0.958613029344 , -0.958612574286 , -0.958612461896 , -0.958612433968 , -0.958612427006 , -0.958612425271 , -0.958612424832 , -0.958612424732]
#  RK4_noadjust_final_error_list=[ -0.958612424733 , -0.958612424693 , -0.95861242469 , -0.958612424691 , -0.95861242469 , -0.958612424692 , -0.958612424688 , -0.958612424696]
#  RK4_adjusted_final_error_list=[ -0.958612424733 , -0.958612424693 , -0.95861242469 , -0.958612424691 , -0.95861242469 , -0.958612424692 , -0.958612424688 , -0.958612424696]

# /r/T/C/T/p/m/test-fig2> cat output.txt | grep VOL3
#  canonical                    =  -4.65783235859
#  RK2_noadjust_final_error_list=[ -4.6578300261 , -4.65783177912 , -4.65783221417 , -4.65783232254 , -4.65783234958 , -4.65783235633 , -4.65783235802 , -4.65783235844]
#  RK2_adjusted_final_error_list=[ -4.6578300261 , -4.65783177912 , -4.65783221417 , -4.65783232254 , -4.65783234958 , -4.65783235633 , -4.65783235802 , -4.65783235844]
#  RK4_noadjust_final_error_list=[ -4.65783235852 , -4.65783235858 , -4.65783235858 , -4.65783235858 , -4.65783235858 , -4.65783235858 , -4.65783235858 , -4.65783235858]
#  RK4_adjusted_final_error_list=[ -4.65783235852 , -4.65783235858 , -4.65783235858 , -4.65783235858 , -4.65783235858 , -4.65783235858 , -4.65783235858 , -4.65783235858]

# /r/T/C/T/p/m/test-fig2> cat output.txt | grep VOL4
#  canonical                    =  -0.748436417459
#  RK2_noadjust_final_error_list=[ -0.74843707901 , -0.748436580988 , -0.74843645808 , -0.748436427549 , -0.748436419941 , -0.74843641804 , -0.74843641757 , -0.748436417443]
#  RK2_adjusted_final_error_list=[ -0.74843707901 , -0.748436580988 , -0.74843645808 , -0.748436427549 , -0.748436419941 , -0.74843641804 , -0.74843641757 , -0.748436417443]
#  RK4_noadjust_final_error_list=[ -0.748436417461 , -0.748436417412 , -0.748436417409 , -0.748436417409 , -0.74843641741 , -0.748436417408 , -0.748436417412 , -0.748436417403]
#  RK4_adjusted_final_error_list=[ -0.748436417461 , -0.748436417412 , -0.748436417409 , -0.748436417409 , -0.74843641741 , -0.748436417408 , -0.748436417412 , -0.748436417403]


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

# time step lengthes: [0.001, 0.0005, 0.00025, 0.000125, 6.25e-05, 3.125e-05, 1.5625e-05] ...

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

print "\n fit RK2 noAdjust:", p2n
print "\n fit RK2 Adjusted:", p2a
print "\n fit RK4 noAdjust:", p4n
print "\n fit RK4 Adjusted:", p4a

plot(ts_list, p2n(ts_list), 'b', label='%s'%p2n)
plot(ts_list, p2a(ts_list), 'b', label='%s'%p2a)
plot(ts_list, p4n(ts_list), 'r', label='%s'%p4n)
plot(ts_list[:5], p4a(ts_list[:5]), 'r', label='%s'%p4a)

legend(loc='upper left')
xlabel("$\log_{10}$(time_step)")
ylabel("$\log_{10}$(error)")
# title("")
savefig('final_error_list.pdf', dpi=600)

