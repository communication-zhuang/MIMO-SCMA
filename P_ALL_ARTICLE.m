%%  
%第三章仿真图
% 1. BsMP与MPA算法在不同信道、不同码本、不同配置参数下的仿真图   ok
% 2. BsMP（1-1）在瑞利信道下，使用华为码本，多接收天线不编码的仿真图  ok
% 3. BsMP（1-1）在瑞利信道下，使用华为码本，多发送天线多接收天线的仿真图。 ok

%第五章仿真图
%1. BsMP在LDPC编码下的SDD仿真图   ok
%2. BsMP在LDPC编码下的IDD仿真图   no


%%  iter_MPA = 8, iter_BsMP = 8(thre = 2), Nt = 1, Nr = 1, Rayli channel, MPA, BsMP, huawei
clc;
clear

SNR = 0 : 2 : 18;
MPA = [0.211790393013100	0.169894366197183	0.140476190476190	0.0775316455696203	0.0479483633010604	0.0258077226162333	0.00916299096828397	0.00299406853300087	0.00122549019607843	0.000437215809723680];
BsMP3 = [0.211790393013100	0.173886883273165	0.139784946236559	0.0837095560571859	0.0503099666189795	0.0296660703637448	0.0107417038137692	0.00453560790685080	0.00196893063583815	0.000978309841024651];
BsMP2 = [0.223467862481315	0.175559947299078	0.162783751493429	0.0968406593406593	0.0709186840471757	0.0427897574123989	0.0237481590574374	0.0118912440070654	0.00683621933621934	0.00334542523176139];
BsMP1 = [0.264492753623188	0.225975975975976	0.243939393939394	0.198412698412698	0.174586776859504	0.156443618339529	0.121026339691190	0.0952119883040936	0.0744485294117647	0.0436507936507937];
s0 = semilogy(SNR, BsMP3, '-d', 'linewidth', 1.5);
set(s0,'Color',[128 41 41]/256)

hold on;

s1 = semilogy(SNR, BsMP2, '-o', 'linewidth', 1.5);
set(s1,'Color',[0 0 256]/256)


s2 = semilogy(SNR, BsMP1, '-<', 'linewidth', 1.5);
set(s2,'Color',[255 62 150]/256)

s3 = semilogy(SNR, MPA, '-h', 'linewidth', 1.5);
set(s3,'Color',[128,0,128]/256)



grid on;
%title('Comparision of FER for Max-Log ','interpreter','latex');
legend('BsMP, $I=8,T=2,d_s=3$','BsMP, $I=8,T=2,d_s=2$','BsMP, $I=8,T=2,d_s=1$','Max-Log, $I=8$','interpreter','latex');
axis([0, 18, 1e-04, 0.4]);
xlabel('SNR[dB]', 'interpreter','latex');
ylabel('BER', 'interpreter','latex');

%%  iter_MPA = 8, iter_BsMP = 8(thre = 2), Nt = 1, Nr = 1, AWGN, MPA, BsMP, huawei
clc;
clear

SNR = 0 : 2 : 12;
MPA = [0.188888888888889	0.117268041237113	0.0671589310829817	0.0272152994115611	0.00709864082895976	0.00180237715096683	0.000351674245230233];
BsMP3 = [0.185873605947955	0.121573828470380	0.0675451092117759	0.0278505325650428	0.00729920237432758	0.00186048413477265	0.000353786954231105];
BsMP2 = [0.186419753086420	0.133045148895293	0.0845388788426763	0.0337203847056064	0.0118772173376523	0.00393597630439736	0.000780118132174301];
BsMP1 = [0.215530303030303	0.199353448275862	0.179218967921897	0.133164414414414	0.0865026099925429	0.0461783439490446	0.0146052212671124];

s0 = semilogy(SNR, BsMP3, '-d', 'linewidth', 1.5);
set(s0,'Color',[128 41 41]/256)

hold on;

s1 = semilogy(SNR, BsMP2, '-o', 'linewidth', 1.5);
set(s1,'Color',[0 0 256]/256)


s2 = semilogy(SNR, BsMP1, '-<', 'linewidth', 1.5);
set(s2,'Color',[255 62 150]/256)

s3 = semilogy(SNR, MPA, '-h', 'linewidth', 1.5);
set(s3,'Color',[128,0,128]/256)



grid on;
%title('Comparision of FER for Max-Log ','interpreter','latex');
legend('BsMP, $I=8,T=2,d_s=3$','BsMP, $I=8,T=2,d_s=2$','BsMP, $I=8,T=2,d_s=1$','Max-Log, $I=8$','interpreter','latex');
axis([0, 12, 1e-04, 0.4]);
xlabel('SNR[dB]', 'interpreter','latex');
ylabel('BER', 'interpreter','latex');

%%  iter_EP = 8, iter_MPA = 8, iter_BsMP = 8(thre = 2), Nt = 1, Nr = 2, Rayli channel, MPA, BsMP, huawei
clc;
clear

SNR = 0 : 2 : 12;
MPA = [0.0779376498800959	0.0401484480431849	0.0135587251276633	0.00482456140350877	0.00105786387968392	0.000224016899453589	5.78915649809806e-05];
BsMP1 = [0.109887005649718	0.0811800172265289	0.0489790286975717	0.0230637060533968	0.00963756177924218	0.00271978829755954	0.000874692237916289];
EP = [0.0765306122448980	0.0437201907790143	0.0207070707070707	0.0105541681531692	0.00643500643500644	0.00498673604734211	0.00388426476416964];
s0 = semilogy(SNR, BsMP1, '-d', 'linewidth', 1.5);
set(s0,'Color',[128 41 41]/256)
hold on;

s3 = semilogy(SNR, MPA, '-h', 'linewidth', 1.5);
set(s3,'Color',[128,0,128]/256)

s2 = semilogy(SNR, EP, '-<', 'linewidth', 1.5);
set(s2,'Color',[255 62 150]/256)

grid on;
%title('Comparision of FER for Max-Log ','interpreter','latex');
legend('BsMP, $I=8,T=2,d_s=1$','Max-Log, $I=8$', 'EP, $I=8$','interpreter','latex');
axis([0, 12, 1e-05, 0.2]);
xlabel('SNR[dB]', 'interpreter','latex');
ylabel('BER', 'interpreter','latex');

%%  iter_EP = 8, iter_MPA = 8, iter_BsMP = 8(thre = 2), Nt = 1, Nr = 4, Rayli channel, MPA, BsMP, EP, uncoded system, huawei
clc;
clear

SNR = -6 : 2 : 6;
MPA = [0.115313653136531 0.0701407211961302	0.0291241496598639 0.00870420829548895	0.00218609865470852	0.000295351661698134	4.21227034351065e-05];
BsMP1 = [0.123517786561265 0.0813008130081301	0.0351465474416294 0.0151353276353276	0.00325605274535236	0.000504983388704319	6.70130727040290e-05];
EP = [0.107716049382716 0.0662280701754386	0.0267489711934156	0.0107049608355091	0.00319846658352963	0.00124916429149513	0.000741139726759528];
s0 = semilogy(SNR, BsMP1, '-d', 'linewidth', 1.5);
set(s0,'Color',[128 41 41]/256)
hold on;

s3 = semilogy(SNR, MPA, '-h', 'linewidth', 1.5);
set(s3,'Color',[128,0,128]/256)

s2 = semilogy(SNR, EP, '-<', 'linewidth', 1.5);
set(s2,'Color',[255 62 150]/256)

grid on;
%title('Comparision of FER for Max-Log ','interpreter','latex');
legend('BsMP, $I=8,T=2,d_s=1$','Max-Log, $I=8$', 'EP, $I = 8$', 'interpreter','latex');
axis([-6, 6, 1e-05, 0.2]);
xlabel('SNR[dB]', 'interpreter','latex');
ylabel('BER', 'interpreter','latex');



%% iter_EP = 6, iter_MPA = 6, iter_BsMP = 6(thre = 2), Nt = 1, Nr = 4, Rayli channel, MPA, BsMP, EP, LDPC coded system, R = 3/4, iter_LDPC = 20, SDD, huawei
%  N0 of coded system is equals to the N0 of uncoded system
clc;
clear

SNR = -7 : 4;
EP4_uncoded = [0.134351851851852	0.108834876543210	0.0843569958847737	0.0617763051722747	0.0428252759509250	0.0294087153475612	0.0192889436327255	0.0109863898659678...
                0.00642928076306577	0.00334581970384440	0.00185802046317890	0.00113840790727953];

EP4_coded = [0.133443072702332	0.104602194787380	0.0752057613168724	0.0477890744775276	0.0245054677384203	0.0118565306621580	0.00543820238859315	0.00165720241242691...
                0.000618703343786150	0.000226669636820528	0.000121711156171420	9.68744429719529e-05];


MPA4_uncoded = [0.141630658436214	0.116705246913580	0.0908539094650206	0.0681995884773663	0.0473968225704337	0.0316083957363557	0.0187613080034057	0.00992895411134405...
                0.00491823190059260	0.00221707224930859	0.000921920795200122 0.000348621420580357];

MPA4_coded = [0.140510973936900	0.113539094650206	0.0833058984910837	0.0554561042524006	0.0296672206394429	0.0150441880860824	0.00496221796509153	0.00150350464815874...
              0.000345552246860083	8.22490868252565e-05	1.67725137035595e-05 3.47887110416672e-06];

BsMP4_uncoded = [0.147384259259259	0.122584876543210	0.0994444444444445	0.0763837448559671	0.0524176954732510	0.0370848355730362	0.0237214289255789	0.0130870705042392...
                 0.00670110644299435	0.00314366409634791	0.00130514098160461 0.000500380200329326];

BsMP4_coded = [0.146409465020576	0.119142661179698	0.0941289437585734	0.0662311385459534	0.0366392318244170	0.0204183206477536	0.00865189151668402	0.00271456195150974...
                0.000647245976646173	0.000134825725370891	2.68609630795029e-05 5.23694357518869e-06];




s6 = semilogy(SNR, EP4_uncoded, '-h', 'linewidth', 1.2);
set(s6,'Color',[128,0,128]/256)
hold on

s2 = semilogy(SNR, MPA4_uncoded, '-s', 'linewidth', 1.2);
set(s2,'Color',[128,0,128]/256)


s9 = semilogy(SNR, BsMP4_uncoded, '-*', 'linewidth', 1.2);
set(s9,'Color',[128,0,128]/256)

s0 = semilogy(SNR, EP4_coded, '-h', 'linewidth', 1.2);
set(s0,'Color',[252 170 103]/256)

s3 = semilogy(SNR, MPA4_coded, '-s', 'linewidth', 1.2);
set(s3,'Color',[252 170 103]/256)


s1 = semilogy(SNR, BsMP4_coded, '-*', 'linewidth', 1.2);
set(s1,'Color',[252 170 103]/256)


grid on;
%title('Comparision of BER for different antennas','interpreter','latex');
legend('GA,$N_t=1,N_r=4$','MPA,$N_t=1,N_r=4$','BsMP,$N_t=1,N_r=4$', 'EP4', ...
        'MPA with LDPC,$N_t=1,N_r=4$','BsMP with LDPC,$N_t=1,N_r=4$','interpreter','latex');
axis([-7, 4, 3e-06, 0.2]);
xlabel('SNR[dB]', 'interpreter','latex');
ylabel('BER', 'interpreter','latex');


%% (1-6-20)iter_MPA = 6, Nt = 1, Nr = 4, Rayli channel, MPA, LDPC coded system, R = 3/4, iter_LDPC = 20, IDD, huawei， dynamic weights
%  N0 of coded system is equals to the N0 of uncoded system
clc;
clear

SNR = -7:2:1;
MPA4_uncoded = [0.141630658436214 0.0908539094650206 0.0473968225704337 0.0187613080034057 0.00491823190059260];
MPA4_SDD = [0.140510973936900	0.0833058984910837 0.0296672206394429 0.00496221796509153	0.000345552246860083];
MPA4_IDD_2_3_10 = [0.138676268861454	0.0761454046639232	0.0211676954732510	0.00290111718165689	0.000245934869715447];
MPA4_IDD_3_2_7 = [0.138744855967078	0.0778086419753086	0.0220226960967702	0.00291031772513254	0.000282207361682670];



s6 = semilogy(SNR, MPA4_uncoded, '-h', 'linewidth', 1.2);
set(s6,'Color',[128,0,128]/256)
hold on

s2 = semilogy(SNR, MPA4_SDD, '-s', 'linewidth', 1.2);
set(s2,'Color',[128,0,128]/256)


s9 = semilogy(SNR, MPA4_IDD_2_3_10, '-*', 'linewidth', 1.2);
set(s9,'Color',[128,0,128]/256)

s0 = semilogy(SNR, MPA4_IDD_3_2_7, '-h', 'linewidth', 1.2);
set(s0,'Color',[252 170 103]/256)



grid on;
%title('Comparision of BER for different antennas','interpreter','latex');
legend('MPA', 'MPA-SDD', 'MPA-IDD1','MPA-IDD2', 'interpreter','latex');
axis([-7, 1, 1e-4, 0.2]);
xlabel('SNR[dB]', 'interpreter','latex');
ylabel('BER', 'interpreter','latex');

%% (1-6-20) iter_MPA = 6, Nt = 1, Nr = 4, AWGN channel, MPA, LDPC coded system, R = 3/4, iter_LDPC = 20, IDD, huawei， dynamic weights
%  N0 of coded system is equals to the N0 of uncoded system
clc;
clear

SNR = -2:5;
MPA4_uncoded = [0.0895550411522634 0.0703463203463204	0.0533892727871010	0.0434442548350398	0.0308033514046328	0.0251197809719370	0.0194316657082615	0.0119749550763702];
MPA4_SDD = [0.107887517146776 0.0903017832647462	0.0725445816186557	0.0520130315500686	0.0276953413527488	0.00746944405912363	0.00170774965013648	0.000221074790140711];
MPA4_IDD = [0.130617283950617 0.117177640603567	0.0906138545953361	0.0481543833395685	0.00886967441297507	0.000782401056749479	3.88254863203828e-05 0];
MPA4_dy_IDD = [0.0818381344307270 0.0779149519890261	0.0719513031550069	0.0371295751954188	0.00426383173296754	0.000279390989138393 1.84110351013070e-05 0];



s6 = semilogy(SNR, MPA4_uncoded, '-h', 'linewidth', 1.2);
set(s6,'Color',[128,0,128]/256)
hold on

s2 = semilogy(SNR, MPA4_SDD, '-s', 'linewidth', 1.2);
set(s2,'Color',[128,0,128]/256)


s9 = semilogy(SNR, MPA4_IDD, '-*', 'linewidth', 1.2);
set(s9,'Color',[128,0,128]/256)

s0 = semilogy(SNR, MPA4_dy_IDD, '-h', 'linewidth', 1.2);
set(s0,'Color',[252 170 103]/256)



grid on;
%title('Comparision of BER for different antennas','interpreter','latex');
legend('MPA', 'MPA-SDD', 'MPA-IDD', 'MPA-dy-IDD','interpreter','latex');
axis([-1, 5, 1e-5, 0.2]);
xlabel('SNR[dB]', 'interpreter','latex');
ylabel('BER', 'interpreter','latex');

































%% iter_EP = 8, iter_MPA = 8, iter_BsMP = 8(thre = 2), Nt = 1, Nr = 4, Rayli channel, MPA, BsMP, EP, LDPC coded system, R = 3/4, iter_LDPC = 20, SDD, huawei
%  N0 of coded system is not equals to the N0 of uncoded system
clc;
clear

SNR = -6 : 2 : 4;
MPA_uncode = [0.115313653136531 0.0701407211961302	0.0291241496598639 0.00870420829548895	0.00218609865470852	0.000295351661698134];
BsMP1_uncode = [0.123517786561265 0.0813008130081301	0.0351465474416294 0.0151353276353276	0.00325605274535236	0.000504983388704319];
EP_uncode = [0.107716049382716 0.0662280701754386	0.0267489711934156	0.0107049608355091	0.00319846658352963	0.00124916429149513	];

MPA_code = [0.154243827160494	0.0926041666666667	0.0303497942386831	0.00502188364507205	0.000327492019903300	1.12842635485263e-05];
BsMP1_code = [0.160277777777778	0.104405864197531	0.0496310047671434	0.00944989106753813	0.000600291928510942	2.74566054471807e-05];
EP_code = [0.149467592592593	0.0827006172839506	0.0297428464289835	0.00471835075493612	0.000328386996514254	6.60891103801560e-05];

s6 = semilogy(SNR, EP_uncode, '-h', 'linewidth', 1.2);
set(s6,'Color',[128,0,128]/256)
hold on

s2 = semilogy(SNR, BsMP1_uncode, '-s', 'linewidth', 1.2);
set(s2,'Color',[128,0,128]/256)


s9 = semilogy(SNR, MPA_uncode, '-*', 'linewidth', 1.2);
set(s9,'Color',[128,0,128]/256)

s0 = semilogy(SNR, EP_code, '-h', 'linewidth', 1.2);
set(s0,'Color',[252 170 103]/256)

s3 = semilogy(SNR, BsMP1_code, '-s', 'linewidth', 1.2);
set(s3,'Color',[252 170 103]/256)


s1 = semilogy(SNR, MPA_code, '-*', 'linewidth', 1.2);
set(s1,'Color',[252 170 103]/256)


grid on;
%title('Comparision of BER for different antennas','interpreter','latex');
legend('EP,$N_t=1,N_r=4$','MPA,$N_t=1,N_r=4$','BsMP,$N_t=1,N_r=4$', 'EP with LDPC,$N_t=1,N_r=4$', ...
        'BsMP with LDPC,$N_t=1,N_r=4$','MPA with LDPC,$N_t=1,N_r=4$','interpreter','latex');
axis([-6, 4, 3e-06, 0.2]);
xlabel('SNR[dB]', 'interpreter','latex');
ylabel('BER', 'interpreter','latex');



%% iter_EP = 12, Nt = 1, Nr = 4, Rayli channel, EP, LDPC coded system, R = 3 / 4, iter_LDPC = 12, IDD, huawei
%  N0 of coded system is equals to the N0 of uncoded system
clc;
clear

SNR = -8 : 1;

EP_1_12_12 = [0.159955418381344	0.131124828532236	0.100216049382716	0.0692249657064472	0.0410528120713306	0.0198320840681952	0.0108848627367146	0.00349598275524201	0.00125783136396089	0.000632937735937834];
EP_2_6_6 = [0.160174897119342	0.131985596707819	0.103264746227709	0.0727743484224966	0.0448148148148148	0.0222691827630099	0.0120587787012382	0.00388307054453538	0.00145354694469026	0.000717828415867072];
EP_3_4_4 = [0.160833333333333	0.133048696844993	0.105740740740741	0.0760425240054870	0.0467729766803841	0.0239158911998418	0.0127067833101682	0.00411869033268210	0.00177693522513380	0.000676575219886946];
EP_4_3_3 = [0.161292866941015	0.134348422496571	0.107397119341564	0.0773593964334705	0.0475274348422497	0.0237355840065031	0.0123772971816050	0.00392355193452587	0.00137300600263563	0.000518424484483147];
EP_6_2_2 = [0.161646090534979	0.135130315500686	0.107870370370370	0.0782853223593964	0.0487722908093278	0.0260239075053890	0.0133437749523985	0.00508368804347186	0.00179185369411894	0.000696511060693859];
EP_12_1_1 = [0.162280521262003	0.136601508916324	0.109173525377229	0.0786248285322359	0.0479766803840878	0.0231591598394181	0.0117799605498984	0.00369974675530231	0.00134010969920085	0.000448084482090461];
s6 = semilogy(SNR, EP_1_12_12, '-h', 'linewidth', 1.2);
set(s6,'Color',[128,0,128]/256)
hold on

s2 = semilogy(SNR, EP_2_6_6, '-s', 'linewidth', 1.2);
set(s2,'Color',[128,0,128]/256)


s9 = semilogy(SNR, EP_3_4_4, '-*', 'linewidth', 1.2);
set(s9,'Color',[128,0,128]/256)

s0 = semilogy(SNR, EP_4_3_3, '-^', 'linewidth', 1.2);
set(s0,'Color',[128,0,128]/256)

s1 = semilogy(SNR, EP_6_2_2, '-o', 'linewidth', 1.2);
set(s1,'Color',[128,0,128]/256)

s5 = semilogy(SNR, EP_12_1_1, '--', 'linewidth', 1.2);
set(s5,'Color',[128,0,128]/256)



grid on;
%title('Comparision of BER for different antennas','interpreter','latex');
legend('EP11212,$N_t=1,N_r=4$','EP266,$N_t=1,N_r=4$','EP344,$N_t=1,N_r=4$', 'EP433,$N_t=1,N_r=4$', ...
    'EP622,$N_t=1,N_r=4$','EP1211,$N_t=1,N_r=4$','interpreter','latex');
axis([-8, 1, 3e-4, 0.2]);
xlabel('SNR[dB]', 'interpreter','latex');
ylabel('BER', 'interpreter','latex');


