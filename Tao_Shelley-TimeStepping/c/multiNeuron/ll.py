canonical1                    =  0.175578230761093
RK2_noadjust_final_error_list1=[ 0.905871657082677 , 0.982679406940519 , 0.048688959563439 , 0.108893854748690 , 0.145775348112375 , 0.161937648514800 , 0.169672260781170 , 0.172953565766778]
RK2_adjusted_final_error_list1=[ 0.156027487621844 , 0.170849404180191 , 0.174410659377428 , 0.175288424862467 , 0.175505144936786 , 0.175560056287306 , 0.175573709228726 , 0.175577103729076]
RK4_noadjust_final_error_list1=[ 0.906503285509976 , 0.982822531786080 , 0.071995251519531 , 0.108903143689566 , 0.145777722981026 , 0.161938246801384 , 0.169672410923413 , 0.172953603358968]
RK4_adjusted_final_error_list1=[ 0.175573571998176 , 0.175577947738283 , 0.175578214130311 , 0.175578230531152 , 0.175578231542195 , 0.175578231622087 , 0.175578231594010 , 0.175578231659159]

canonical2                    =  0.175578230760999
RK2_noadjust_final_error_list2=[ 0.905871657082677 , 0.982679406940519 , 0.048688959563439 , 0.108893854748690 , 0.145775348112375 , 0.161937648514800 , 0.169672260781170 , 0.172953565766778]
RK2_adjusted_final_error_list2=[ 0.156027487621844 , 0.170849404180182 , 0.174410659377428 , 0.175288424862448 , 0.175505144936757 , 0.175560056287306 , 0.175573709228735 , 0.175577103729019]
RK4_noadjust_final_error_list2=[ 0.906503285509976 , 0.982822531786080 , 0.071995251519531 , 0.108903143689566 , 0.145777722981026 , 0.161938246801384 , 0.169672410923413 , 0.172953603358968]
RK4_adjusted_final_error_list2=[ 0.175573571998176 , 0.175577947738273 , 0.175578214130311 , 0.175578230531142 , 0.175578231542195 , 0.175578231622115 , 0.175578231594039 , 0.175578231659159]

canonical3                    =  0.175578230761093
RK2_noadjust_final_error_list3=[ 0.905871657082677 , 0.982679406940519 , 0.048688959563439 , 0.108893854748690 , 0.145775348112375 , 0.161937648514800 , 0.169672260781170 , 0.172953565766778]
RK2_adjusted_final_error_list3=[ 0.156027487621844 , 0.170849404180182 , 0.174410659377419 , 0.175288424862448 , 0.175505144936757 , 0.175560056287306 , 0.175573709228726 , 0.175577103729085]
RK4_noadjust_final_error_list3=[ 0.906503285509976 , 0.982822531786080 , 0.071995251519531 , 0.108903143689566 , 0.145777722981026 , 0.161938246801384 , 0.169672410923413 , 0.172953603358968]
RK4_adjusted_final_error_list3=[ 0.175573571998147 , 0.175577947738254 , 0.175578214130273 , 0.175578230531142 , 0.175578231542195 , 0.175578231622096 , 0.175578231594010 , 0.175578231659178]

canonical4                    =  0.175578230760999
RK2_noadjust_final_error_list4=[ 0.905871657082677 , 0.982679406940519 , 0.048688959563439 , 0.108893854748690 , 0.145775348112375 , 0.161937648514800 , 0.169672260781170 , 0.172953565766778]
RK2_adjusted_final_error_list4=[ 0.156027487621844 , 0.170849404180182 , 0.174410659377428 , 0.175288424862448 , 0.175505144936757 , 0.175560056287306 , 0.175573709228678 , 0.175577103729066]
RK4_noadjust_final_error_list4=[ 0.906503285509976 , 0.982822531786080 , 0.071995251519531 , 0.108903143689566 , 0.145777722981026 , 0.161938246801384 , 0.169672410923413 , 0.172953603358968]
RK4_adjusted_final_error_list4=[ 0.175573571998147 , 0.175577947738273 , 0.175578214130301 , 0.175578230531142 , 0.175578231542195 , 0.175578231622087 , 0.175578231594010 , 0.175578231659159]


for i in range(len(RK2_noadjust_final_error_list1)):
    RK2_noadjust_final_error_list1[i] = abs(RK2_noadjust_final_error_list1[i] - canonical1)
    RK2_noadjust_final_error_list2[i] = abs(RK2_noadjust_final_error_list2[i] - canonical2)
    RK2_noadjust_final_error_list3[i] = abs(RK2_noadjust_final_error_list3[i] - canonical3)
    RK2_noadjust_final_error_list4[i] = abs(RK2_noadjust_final_error_list4[i] - canonical4)
    if (RK2_noadjust_final_error_list1[i]>0.5): RK2_noadjust_final_error_list1[i] = 1 - RK2_noadjust_final_error_list1[i]
    if (RK2_noadjust_final_error_list2[i]>0.5): RK2_noadjust_final_error_list2[i] = 1 - RK2_noadjust_final_error_list2[i]
    if (RK2_noadjust_final_error_list3[i]>0.5): RK2_noadjust_final_error_list3[i] = 1 - RK2_noadjust_final_error_list3[i]
    if (RK2_noadjust_final_error_list4[i]>0.5): RK2_noadjust_final_error_list4[i] = 1 - RK2_noadjust_final_error_list4[i]

for i in range(len(RK2_adjusted_final_error_list1)):
    RK2_adjusted_final_error_list1[i] = abs(RK2_adjusted_final_error_list1[i] - canonical1)
    RK2_adjusted_final_error_list2[i] = abs(RK2_adjusted_final_error_list2[i] - canonical2)
    RK2_adjusted_final_error_list3[i] = abs(RK2_adjusted_final_error_list3[i] - canonical3)
    RK2_adjusted_final_error_list4[i] = abs(RK2_adjusted_final_error_list4[i] - canonical4)
    if (RK2_adjusted_final_error_list1[i]>0.5): RK2_adjusted_final_error_list1[i] = 1 - RK2_adjusted_final_error_list1[i]
    if (RK2_adjusted_final_error_list2[i]>0.5): RK2_adjusted_final_error_list2[i] = 1 - RK2_adjusted_final_error_list2[i]
    if (RK2_adjusted_final_error_list3[i]>0.5): RK2_adjusted_final_error_list3[i] = 1 - RK2_adjusted_final_error_list3[i]
    if (RK2_adjusted_final_error_list4[i]>0.5): RK2_adjusted_final_error_list4[i] = 1 - RK2_adjusted_final_error_list4[i]

for i in range(len(RK4_noadjust_final_error_list1)):
    RK4_noadjust_final_error_list1[i] = abs(RK4_noadjust_final_error_list1[i] - canonical1)
    RK4_noadjust_final_error_list2[i] = abs(RK4_noadjust_final_error_list2[i] - canonical2)
    RK4_noadjust_final_error_list3[i] = abs(RK4_noadjust_final_error_list3[i] - canonical3)
    RK4_noadjust_final_error_list4[i] = abs(RK4_noadjust_final_error_list4[i] - canonical4)
    if (RK4_noadjust_final_error_list1[i]>0.5): RK4_noadjust_final_error_list1[i] = 1 - RK4_noadjust_final_error_list1[i]
    if (RK4_noadjust_final_error_list2[i]>0.5): RK4_noadjust_final_error_list2[i] = 1 - RK4_noadjust_final_error_list2[i]
    if (RK4_noadjust_final_error_list3[i]>0.5): RK4_noadjust_final_error_list3[i] = 1 - RK4_noadjust_final_error_list3[i]
    if (RK4_noadjust_final_error_list4[i]>0.5): RK4_noadjust_final_error_list4[i] = 1 - RK4_noadjust_final_error_list4[i]

for i in range(len(RK4_adjusted_final_error_list1)):
    RK4_adjusted_final_error_list1[i] = abs(RK4_adjusted_final_error_list1[i] - canonical1)
    RK4_adjusted_final_error_list2[i] = abs(RK4_adjusted_final_error_list2[i] - canonical2)
    RK4_adjusted_final_error_list3[i] = abs(RK4_adjusted_final_error_list3[i] - canonical3)
    RK4_adjusted_final_error_list4[i] = abs(RK4_adjusted_final_error_list4[i] - canonical4)
    if (RK4_adjusted_final_error_list1[i]>0.5): RK4_adjusted_final_error_list1[i] = 1 - RK4_adjusted_final_error_list1[i]
    if (RK4_adjusted_final_error_list2[i]>0.5): RK4_adjusted_final_error_list2[i] = 1 - RK4_adjusted_final_error_list2[i]
    if (RK4_adjusted_final_error_list3[i]>0.5): RK4_adjusted_final_error_list3[i] = 1 - RK4_adjusted_final_error_list3[i]
    if (RK4_adjusted_final_error_list4[i]>0.5): RK4_adjusted_final_error_list4[i] = 1 - RK4_adjusted_final_error_list4[i]


RK2_noadjust_final_error_list=[]
RK2_adjusted_final_error_list=[]
RK4_noadjust_final_error_list=[]
RK4_adjusted_final_error_list=[]

for i in range(len(RK4_adjusted_final_error_list1)):
    RK2_noadjust_final_error_list.append( (RK2_noadjust_final_error_list1[i] + RK2_noadjust_final_error_list2[i] + RK2_noadjust_final_error_list3[i] + RK2_noadjust_final_error_list4[i])/4.0 )
    RK2_adjusted_final_error_list.append( (RK2_adjusted_final_error_list1[i] + RK2_adjusted_final_error_list2[i] + RK2_adjusted_final_error_list3[i] + RK2_adjusted_final_error_list4[i])/4.0 )
    RK4_noadjust_final_error_list.append( (RK4_noadjust_final_error_list1[i] + RK4_noadjust_final_error_list2[i] + RK4_noadjust_final_error_list3[i] + RK4_noadjust_final_error_list4[i])/4.0 )
    RK4_adjusted_final_error_list.append( (RK4_adjusted_final_error_list1[i] + RK4_adjusted_final_error_list2[i] + RK4_adjusted_final_error_list3[i] + RK4_adjusted_final_error_list4[i])/4.0 )


ts_list=[]
for i in range(8): ts_list.append(0.001 / 2**i)

#  figure()
#  plot(ts_list, RK2_noadjust_final_error_list, '+b', label='RK2 Org')
#  plot(ts_list, RK2_adjusted_final_error_list, '.b', label='RK2 Adj')
#  plot(ts_list, RK4_noadjust_final_error_list, 'xr', label='RK4 Org')
#  plot(ts_list, RK4_adjusted_final_error_list, '*r', label='RK4 Adj')

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

