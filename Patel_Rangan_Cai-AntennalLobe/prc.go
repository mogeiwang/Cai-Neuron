package main

// reproduction of Patel's work
// awesome mw, Aug 20, 2015, init, NYU , AbuDhabi
//           , Sep 29, 2015, Ver1, SJTU, Shanghai
//           , Dec 01, 2015, Ver2, SJTU, Shanghai
//           , Dec 20, 2015, Ver3, SJTU, Shanghai
// brilliant , Dec 23, 2015, init, SJTU, Shanghai
//           , Dec 26, 2015, Ver1, SJTU, Shanghai
//           , Dec 30, 2015, Ver2, SJTU, Shanghai
//           , Jan 04, 2016, Ver3, SJTU, Shanghai
// A 1.0.1 update:
// a) init_matrix has its own rand.Seed now
// b) xN_V_list, xN_stim_list sampling rate reset to 1item per ms
// c) add running watch stop settings which can be set in yaml file
//        if_running >0 run, <0 pause
// d) if realtime plot or not can be set in yaml file
//        if_realtime_plot >0 run, <0 pause
// A 1.0.2 update:
// a) do not clear and reset config.yaml, if it exists before running!
// b) show config.yaml each time telling time
// A 1.1.0 update:
// a) change xN_V_list, xN_stim_list's dimension to [time_step][neuron_id]
// b) change if_xxx from float64 to int64
// A 1.2.0 update:
// a) try to change some parameters
// A 1.2.1 update:
// a) can set the matrix_seed in command-line now
// A 1.2.2 update:
// a) change xN_number and stimulus_xN_num
// A 1.3 update:
// a) add scale for changing GABA, nACH and slow current strength
//        *scale* vars and scale_synapse_currents() func
// b) add remove coupling things
//        if_*_coupled vars and remove_coupling_matrix() func
//        >0: true; <0: negative; =0 uncertain
// c) add var use_exist_config to control read from or rewrite config.yaml
// d) add usage() function
// A 1.4 update:
// a) plot fluctions of two given neuron ids (one is for PN, and another LN)
// A 1.4.1 update:
// a) add function proc_filenames() to add an appendix to plotting temp filenames
//                                  so that many processes can be run together...
// A 1.5 update:
// a) use E_{Na,K,Ca,A,B,GABA,nACH,slow}_xN to replace the previous hard codes
// b) have the values of parameters re-checked...
// A 2.0 update:
// a) add cancel multi period control. We will get only 1 loop in 1 run from now.
// b) all if_... things changed to >0 True, <=0 False
// A 2.1 update:
// a) do not use period things (do not mod%), when period_len >= time_end
// A 2.1.1 update:
// a) informs users about using old or new config.yaml
// b) use stim_decay (It has been ignored previously!)
// A 2.1.2 update:
// a) ORN stimulus set to 0, after stim_decay
// A 2.2 update:
// a) change parameter values
//        r1, r2, r3, r4: use program values, instead of paper values.
// A 2.3 update:
// a) change ri parameter values back to paper values
// b) print all scale and decouple information
// A 2.4 update:
// a) resort some definations.
// A 2.5 update:
// a) decay 4S:
//    0-nothing-1-stimulus begin-1.5-normal stimulus-3.5-stimulus decay-7.5-nothing-10
// A 2.5.1 update:
// a) del period_len, no more than 1 period/loop is used now...
// A 2.5.2 update:
// a) restore period_len
// b) change stimulus rise and decay time.
// c) stimulate 30 PN & 10 LN.
// A 2.6 update:
// a) support dir/folder things
// b) by default: do NOT realtim plot, do NOT plot fluctions...
// A 2.7 update:
// a) rearrange stimulate time:
//    0-nothing-1.5-stim rise-2.0-normal stim-4.0-stim decay-6.0-nothing-10
// b) Default file dirs also changed.
// A 2.8 update:
// a) set scale and if_._couple in config file
// b) remove g_{GABA, nACH, slow}_{LN/PN} from config file
// c) use_exist_config set to true by default,
//        use config.yaml parameter values prior to source code!!!
// A 3.0 update:
// a) support odor_{PN/LN}list
// b) change stim_{start, end, rise, decay}_time to ratios
// c) del temp plot files in dest_plot()
// d) always init/dest plots
// e) remove random-seeds from function proc_filenames()
// f) print more variables
// A 3.0.0 update:
// !) [re]check. for testing different run-time, different stim-time...
//        mind the long trial time (20S) and temporal stim periods...
// a) default matrix seed is the constant 0
// B 0.0 update:
// a) use seperated LN-slow>PN and LN-GABA>PN synapse
// b) print matrixes in the program startup
// B 0.1 update:
// a) use seperated LN-slow>PN and LN-GABA>PN connection probablity .neuron._.synapse._prob
// b) use seperated LN-slow>PN and LN-GABA>PN connection control, if_.neuron._.synapse.ed
// B 0.2 update:
// a) adjusting coupling probablities
// b) write coupling probs in config.
// B 0.3 update:
// a) add slow GABA overlap control: if_slowGABA_separate
// b) print if slow GABA overlapped on screen and config.
// c) add neuron_num_coefs to adjust neuron numbers...
// B 0.4 update:
// a) count each neuron's spike number
// B 1.0 update:
// a) save coupling matrix in files
// B 1.1 update:
// a) save spiking rates to files
// B 1.2 update:
// a) print the number of presynaptic neurons for each postsynaptic neuron
// B 1.3 update:
// a) readin coupling matrixes
// B 1.4 update:
// a) plot synapse currents
// B 1.5 update:
// a) display matrixes read from csv.
// B 1.6 update:
// a) plot stim fluction...
// B 1.7 update:
// a) save synapse current to file...
// B 1.8 update:
// a) change some var name.
// B 2.0 update:
// a) bug fix, and more...
// B 2.1 update:
// a) if plot PN /LN things or not???
// B 2.2 update:
// a) plot synapse rasters
// B 2.3 update:
// a) bug fix, ...
// B 2.4 update:
// a) add supports for noise-free / posson-random stimulus
// B 3.0 update:
// a) fix some mapset related codes to make it run faster...
// b) introduce if_init_rt_plots to control gnuplot initialization for later potential use
// c) use if_runtime_trace to control the trace of running of codes or not
// B 3.0.1 update:
// a) plot_fluct_xN_id changed to consts
// b) vars ..._Vfluct_... changed to ..._vFluct_...
// c) vars ..._Vraster... changed to ..._vRaster...
// B 3.1 update:
// a) do NOT print spike number of presynaptic neuron & postsynapric neuron any more
// b) thresholds etc. defined as consts
// B 3.1.1 update:
// a) introduce var freq_scoping to control the calculation range
// b) clean period_len - time_end related codes in get_stim_...()
// B 3.1.2 update:
// a) auto-update coupling rates when neuron_num_coefs is changed.

import (
	"encoding/csv"
	"fmt"
	"github.com/deckarep/golang-set"
	"github.com/sabhiram/go-tracey"
	"github.com/sbinet/go-gnuplot"
	"gopkg.in/yaml.v2"
	"io/ioutil"
	"math"
	"math/rand"
	"os"
	"runtime/pprof"
	"strconv"
	"strings"
	"time"
)

const (
	neuron_num_coefs = 1 // 1: 90-30-30-10;;; 0.1; 0.5; 1; 5; 10; ...
	if_readin_matrix = 0 // readin matrix or randomly set?
	if_init_rt_plots = 1 // >0 init gnuplot, <=0 do not init.
	if_runtime_trace = 0
	plot_fluct_PN_id = 0 // which PN fluction to plot?
	plot_fluct_LN_id = 0
	freq_scoping     = "stimulated" // noStim, rising, stimulated, decay, noStim
	// ...
	ms_per_second int64 = 1000               // 1000 ms = 1 S
	time_end      int64 = 10 * ms_per_second // run how many ms; always start from 0! !!!
	period_len    int64 = time_end           // set to "time_end" to prevent stimulating periodically
	PN_number     int64 = 90 * neuron_num_coefs
	LN_number     int64 = 30 * neuron_num_coefs
	stim_PN_num   int64 = 30 * neuron_num_coefs // how many PNs will receive stimulus
	stim_LN_num   int64 = 10 * neuron_num_coefs
	// ...
	ORN_number  int64   = 200
	stim_onset  float64 = 0.15 * float64(time_end) // 0.15 stimulus starts at this moment !!! in ms
	stim_rise   float64 = 0.05 * float64(time_end) // 0.05 rise time of stimulus input onset (1.5-2)
	stim_offset float64 = 0.40 * float64(time_end) // 0.40 withdraw stimulus at this moment (time in ms)
	stim_decay  float64 = 0.20 * float64(time_end) // 0.20 decay time (Inf) of stimulus input offset (4.0-6)
	// ...
	use_exist_config bool    = true                        // read the config file [or rewrite it]?
	plot_file_dir    string  = "/mnt/tmpDisk/"             // file_dir
	save_file_dir    string  = "./"                        // "/mnt/tmpDisk/" // file_dir
	fading_rate      float64 = 0.99800                     // exp(-0.01/5)=0.998
	click_per_ms     int64   = 100                         // time steps (click/iteration) per ms
	click_num_total  int64   = click_per_ms * time_end     // time steps (clicks) in total
	time_stepLen     float64 = 1.0 / float64(click_per_ms) // each time step 0.01ms (corresponding to click_per_ms)
	// ...
	const_scale_0   float64 = 0
	const_scale_050 float64 = 0.50
	const_scale_025 float64 = 0.25
	const_scale_033 float64 = 0.33
	const_scale_1   float64 = 1.0
	const_scale_2   float64 = 2.0
	const_scale_3   float64 = 3.0
	const_scale_4   float64 = 4.0
)

var (
	// scale & remove-coupling things.
	scale_slow_PN float64 = const_scale_1 // change this one
	scale_GABA_PN float64 = const_scale_1 // change this one
	scale_GABA_LN float64 = const_scale_1 // NO need change this
	scale_nACH_PN float64 = const_scale_1 // this can be changed
	scale_nACH_LN float64 = const_scale_1 // this can be changed
	// ...
	if_LN2PN_slowed int64 = 1 // if has LN-2-PN links ; slow
	if_LN2PN_GABAed int64 = 1 // if has LN-2-PN links ; GABA
	if_LN2LN_GABAed int64 = 1 // if has LN-2-LN links
	if_PN2PN_nACHed int64 = 1 // if has PN-2-PN links ; >0: true; <0: false
	if_PN2LN_nACHed int64 = 1 // if has PN-2-LN links
	// ...
	LN2PN_slow_prob float64 = 0.15 / neuron_num_coefs
	LN2PN_GABA_prob float64 = 0.15 / neuron_num_coefs
	LN2LN_GABA_prob float64 = 0.25 / neuron_num_coefs
	PN2PN_nACH_prob float64 = 0.10 / neuron_num_coefs
	PN2LN_nACH_prob float64 = 0.10 / neuron_num_coefs
	// if run, if plot...
	if_running          int64 = 1 // >0 set to run , <0 set to pause.
	if_slowGABA_overlap int64 = 1 // LN2PN slow == GABA ??
	if_random_stim      int64 = 1 // possion input? or stable ones?
	if_realtime_plot    int64 = 1 // >0 set to plot, <0 set to pause.
	if_PN_plot          int64 = 1 // >0 set to plot, <0 set to pause.
	if_LN_plot          int64 = 1 // >0 set to plot, <0 set to pause.
	if_stimRaster_plot  int64 = 1 // >0 set to plot, <0 set to pause.
	if_vRaster_plot     int64 = 1 // >0 set to plot, <0 set to pause.
	if_synaRaster_plot  int64 = 1 // >0 set to plot, <0 set to pause.
	if_vFluct_plot      int64 = 1 // >0 set to plot, <0 set to pause. voltages of PNs and LNs
	if_stimFluct_plot   int64 = 1 // >0 set to plot, <0 set to pause. stim-In...
	if_synaFluct_plot   int64 = 1 // >0 set to plot, <0 set to pause. synapse currents
	//   \--- above variables are set in config file.
	BG_input_rate         float64 = 3.50000 // 3.4711...
	ORN_input_rate        float64 = 0.03500 // 0.0357
	BG_input_strength_PN  float64 = 0.06540
	BG_input_strength_LN  float64 = 0.00010
	ORN_input_strength_PN float64 = 0.01743
	ORN_input_strength_LN float64 = 0.01667
	// ...
	LN_stim_baseline float64 = -0.00174 // avg stim when there is only background
	LN_stim_plusORNs float64 = -0.58715 // avg stim when there are also ORN input
	PN_stim_baseline float64 = -1.14423
	PN_stim_plusORNs float64 = -1.79143
	// more configurations:
	matrix_seed int64 = 0 // rand.seed() in init_matrix()
	odor_slip   int64 = 0 // how many stimulated PNs are chosen from right (max IDs)
	odor_PN_set mapset.Set
	odor_LN_set mapset.Set
	tracy_opt   = tracey.Options{DisableTracing: if_runtime_trace <= 0}
	Exit, Enter = tracey.New(&tracy_opt)
)

////////////////////////////////////////////////////////////
// For this part``` ////////////////////////////////////////
///////////////// define the neuron classes ////////////////
////////////////////////////////////////////////////////////

const (
	V_rest_PN      = -75.0
	V_rest_LN      = -60.0
	V_threshold_PN = 0
	V_threshold_LN = -20
	// ...
	C_m_PN float64 = 1.0
	g_L_PN float64 = 0.3
	E_L_PN float64 = 64.0 // E_... things are timed -1 in this implement
	C_m_LN float64 = 1.0
	g_L_LN float64 = 0.3
	E_L_LN float64 = 50.0
	// ...
	g_K_PN  float64 = 3.6
	g_Na_PN float64 = 120.0
	g_A_PN  float64 = 1.43
	E_K_PN  float64 = 87.0
	E_Na_PN float64 = -40.0
	E_A_PN  float64 = 87.0
	// ...
	g_K_LN  float64 = 36.0
	g_Ca_LN float64 = 5.0
	g_B_LN  float64 = 0.045
	E_K_LN  float64 = 95.0
	E_Ca_LN float64 = -140.0
	E_B_LN  float64 = 95.0
)

var (
	g_GABA_PN float64 = 0.36  // LN --[GABA]-> PN
	g_nACH_PN float64 = 0.009 // PN --[nACH]-> PN
	g_slow_PN float64 = 0.36  // LN --[slow]-> PN
	g_GABA_LN float64 = 0.3   // LN --[GABA]-> LN
	g_nACH_LN float64 = 0.045 // PN --[nACH]-> LN
	// ...
	E_GABA_PN float64 = 70.0
	E_nACH_PN float64 = 0.0
	E_slow_PN float64 = 95.0
	E_GABA_LN float64 = 70.0
	E_nACH_LN float64 = 0.0
	// ...
	CLOCK float64 = 0.0
	err   error   = nil
)

type PNtype struct { //PN struct
	// inner states
	m_K       float64 // PN, LN (inner state)
	h_K       float64 // PN, LN (inner state)
	m_Na      float64 // PN (inner state)
	h_Na      float64 // PN (inner state)
	m_A       float64 // PN (inner state)
	h_A       float64 // PN (inner state)
	O_nACH    float64 // PN, nACH synapse, feed to other neurons
	V         float64 // V (memberane potential)
	lastSpike float64 // last spike time
	// coupling currents
	I_K    float64 // inner current
	I_Na   float64 // inner current
	I_A    float64 // inner current
	I_GABA float64 // synapse in synapse current (from LN)
	I_nACH float64 // synapse in synapse current (from PN)
	I_slow float64 // synapse in from LN
	// background and ORN stimulus
	stim float64 // input stimulations
}

type LNtype struct { //LN struct
	// inner states
	m_K       float64 // PN, LN (inner state)
	h_K       float64 // PN, LN (inner state)
	m_Ca      float64 // LN (inner state)
	h_Ca      float64 // LN (inner state)
	m_B       float64 // LN B; caK (inner state)
	h_B       float64 // LN (inner state)
	F_Ca      float64 // LN (inner state, [Ca])
	O_GABA    float64 // LN, GABA synapse, feed to other neurons
	G_slow    float64 // LN, gSlow, feed to other neurons
	R_slow    float64 // LN, rSlow, feed to other neurons
	V         float64 // V (memberane potential)
	lastSpike float64 // last spike
	// coupling currents
	I_K    float64 // inner current
	I_Ca   float64 // inner current
	I_B    float64 // inner current
	I_GABA float64 // synapse in synapse current (from LN)
	I_nACH float64 // synapse in synapse current (from PN)
	// background and ORN stimulus
	stim float64 // input current
}

type neuronFunctions interface { // mainly for collection,
	dfdt_m_K() float64    //PN, LN
	dfdt_h_K() float64    //PN, LN
	dfdt_m_Na() float64   //PN
	dfdt_h_Na() float64   //PN
	dfdt_m_A() float64    //PN
	dfdt_h_A() float64    //PN
	dfdt_m_Ca() float64   //LN
	dfdt_h_Ca() float64   //LN
	dfdt_m_B() float64    //LN
	dfdt_h_B() float64    //LN
	dfdt_F_Ca() float64   //LN
	dfdt_O_nACH() float64 //PN
	dfdt_O_GABA() float64 //LN
	dfdt_G_slow() float64 //LN
	dfdt_R_slow() float64 //LN
	dfdt_V() float64      //PN, LN
}

// PN methods:
func (x PNtype) dfdt_m_K() float64 {
	var alpha, beta, ret float64
	defer Exit(Enter("$FN()"))
	alpha = (0.01*x.V + 0.5) / (1.0 - math.Exp(-0.1*x.V-5))
	beta = 0.125 * math.Exp(-x.V/80.0-0.75)
	ret = alpha - alpha*x.m_K - beta*x.m_K
	return ret
}

func (x PNtype) dfdt_h_K() float64 {
	defer Exit(Enter("$FN()"))
	return 0.0
}

func (x PNtype) dfdt_m_Na() float64 {
	var alpha, beta, ret float64
	defer Exit(Enter("$FN()"))
	alpha = (0.1*x.V + 3.5) / (1 - math.Exp(-0.1*x.V-3.5))
	beta = 4.0 * math.Exp(-x.V/18.0-10.0/3.0)
	ret = alpha*(1-x.m_Na) - beta*x.m_Na
	return ret
}

func (x PNtype) dfdt_h_Na() float64 {
	var alpha, beta, ret float64
	defer Exit(Enter("$FN()"))
	alpha = 0.07 * math.Exp(-0.05*x.V-3.0)
	beta = 1.0 / (math.Exp(-0.1*x.V-3.0) + 1.0)
	ret = alpha*(1-x.h_Na) - beta*x.h_Na
	return ret
}

func (x PNtype) dfdt_m_A() float64 {
	var mInf, tau, ret float64
	defer Exit(Enter("$FN()"))
	mInf = 1.0 / (1.0 + math.Exp(-x.V/8.5-60.0/8.5))
	tau = 0.27/(math.Exp(x.V/19.7+35.8/19.7)+math.Exp(-x.V/12.7-79.7/12.7)) + 0.1
	ret = (mInf - x.m_A) / tau
	return ret
}

func (x PNtype) dfdt_h_A() float64 {
	var hInf, tau, ret float64
	defer Exit(Enter("$FN()"))
	hInf = 1.0 / (1.0 + math.Exp(x.V/6.0+13.0))
	if x.V < -63.0 {
		tau = 0.27 / (math.Exp(0.2*x.V+9.2) + math.Exp(-x.V/37.5-238.0/37.5))
	} else {
		tau = 5.1
	}
	ret = (hInf - x.h_A) / tau
	return ret
}

func (x PNtype) dfdt_O_nACH() float64 {
	var h1, h2, t1, t2, ret float64
	defer Exit(Enter("$FN()"))
	t1 = x.lastSpike + 0.3 - CLOCK
	t2 = CLOCK - x.lastSpike
	h1 = 1.0 / (1.0 + math.Exp(-100*t1))
	h2 = 1.0 / (1.0 + math.Exp(-100*t2))
	ret = 5.0*h1*h2 - 5.0*h1*h2*x.O_nACH - 0.2*x.O_nACH
	return ret
}

func (x PNtype) dfdt_V() float64 {
	var ret float64
	defer Exit(Enter("$FN()"))
	ret = (-g_L_PN*x.V - g_L_PN*E_L_PN - x.I_Na - x.I_K - x.I_A - x.I_GABA - x.I_nACH - x.stim - x.I_slow) / C_m_PN
	return ret
}

// LN methods:
func (x LNtype) dfdt_m_K() float64 {
	var A, B, ret float64
	defer Exit(Enter("$FN()"))
	A = (-0.70 - 0.02*x.V) / (math.Exp(-7.0-0.2*x.V) - 1.0)
	B = 0.5 * math.Exp(-1-0.025*x.V)
	ret = (A - A*x.m_K - B*x.m_K) / 4.65
	return ret
}

func (x LNtype) dfdt_h_K() float64 {
	defer Exit(Enter("$FN()"))
	return 0.0
}

func (x LNtype) dfdt_m_Ca() float64 {
	var mInf, tau, ret float64
	defer Exit(Enter("$FN()"))
	mInf = 1.0 / (1.0 + math.Exp(-x.V/6.5-20.0/6.5))
	tau = 1.42 + 0.014*x.V
	ret = (mInf - x.m_Ca) / tau
	return ret
}

func (x LNtype) dfdt_h_Ca() float64 {
	var hInf, tau, ret float64
	defer Exit(Enter("$FN()"))
	hInf = 1.0 / (1.0 + math.Exp(x.V/12.0+25.0/12.0))
	tau = 0.3*math.Exp(x.V/13.0-40.0/13.0) + 0.002*math.Exp(-x.V/29.0+60.0/29.0)
	ret = (hInf - x.h_Ca) / tau
	return ret
}

func (x LNtype) dfdt_m_B() float64 {
	var mInf, tau, ret float64
	defer Exit(Enter("$FN()"))
	mInf = x.F_Ca / (x.F_Ca + 2.0)
	tau = 100.0 / (x.F_Ca + 2.0)
	ret = (mInf - x.m_B) / tau
	return ret
}

func (x LNtype) dfdt_h_B() float64 {
	defer Exit(Enter("$FN()"))
	return 0.0
}

func (x LNtype) dfdt_F_Ca() float64 {
	var ret float64
	defer Exit(Enter("$FN()"))
	ret = -(2e-4)*x.I_Ca - x.F_Ca/150.0 + 1.6e-6
	return ret
}

func (x LNtype) dfdt_O_GABA() float64 {
	var temp, ret float64
	defer Exit(Enter("$FN()"))
	temp = 1.0 / (1.0 + math.Exp(-x.V/1.5-20.0/1.5))
	ret = 10.0*temp - 10.0*temp*x.O_GABA - 0.16*x.O_GABA
	return ret
}

func (x LNtype) dfdt_G_slow() float64 {
	var ret, r3, r4 float64 = 0, 0.1000, 0.0330 // !!! r4: 0.033(paper); 0.06(program)
	defer Exit(Enter("$FN()"))
	ret = r3*x.R_slow - r4*x.G_slow
	return ret
}

func (x LNtype) dfdt_R_slow() float64 {
	var ret, h1, h2, t1, t2 float64
	var r1, r2 float64 = 0.5000, 0.0013 // !!! r1, r2: 0.5, 0.0013(paper); 1.0, 0.0025(pragram)
	defer Exit(Enter("$FN()"))
	t1 = x.lastSpike + 0.3 - CLOCK
	t2 = CLOCK - x.lastSpike
	h1 = 1.0 / (1.0 + math.Exp(-100*t1))
	h2 = 1.0 / (1.0 + math.Exp(-100*t2))
	ret = 0.5*r1*h1*h2 - 0.5*r1*h1*h2*x.R_slow - r2*x.R_slow
	return ret
}

func (x LNtype) dfdt_V() float64 {
	var ret float64
	defer Exit(Enter("$FN()"))
	ret = (-g_L_LN*x.V - g_L_LN*E_L_LN - x.I_Ca - x.I_K - x.I_B - x.I_GABA - x.I_nACH - x.stim) / C_m_LN
	return ret
}

////////////////////////////////////////////////////////////
// For this part``` ////////////////////////////////////////
///////////////// define the coupled iteration and more ////
////////////////////////////////////////////////////////////

var (
	PN_reactor     [PN_number][2]PNtype             // PNs' last 2 states&currents structs
	LN_reactor     [LN_number][2]LNtype             // LNs' last 2 states&currents structs
	PN_spike_freq  [PN_number]float64               // how many spike for that PN per second
	LN_spike_freq  [LN_number]float64               // how many spike for that LN per second
	PN_V_list      [time_end + 1][PN_number]float64 // PNs' all historical V list
	LN_V_list      [time_end + 1][LN_number]float64 // LNs' all historical V list
	PN_stim_list   [time_end + 1][PN_number]float64 // PNs' all historical stim list
	LN_stim_list   [time_end + 1][LN_number]float64 // LNs' all historical stim list
	PN_slow_list   [time_end + 1][PN_number]float64 // PNs' all historical feedin slow currents
	PN_GABA_list   [time_end + 1][PN_number]float64 // PNs' all historical feedin GABA currents
	LN_GABA_list   [time_end + 1][LN_number]float64 // LNs' all historical feedin GABA currents
	PN_nACH_list   [time_end + 1][PN_number]float64 // PNs' all historical feedin nACH currents
	LN_nACH_list   [time_end + 1][LN_number]float64 // LNs' all historical feedin nACH currents
	LN2PN_slow_mat [LN_number][PN_number]float64    // coupling matrix
	LN2PN_GABA_mat [LN_number][PN_number]float64    // coupling matrix
	LN2LN_GABA_mat [LN_number][LN_number]float64    // coupling matrix
	PN2PN_nACH_mat [PN_number][PN_number]float64    // coupling matrix
	PN2LN_nACH_mat [PN_number][LN_number]float64    // coupling matrix
)

func readin_matrix() {
	var i, j int64
	defer Exit(Enter("$FN()"))
	if if_readin_matrix <= 0 {
		return
	}
	// ...
	if exist_file(LN2PN_slow_csvName) {
		LN2PN_slow_csv, err = os.Open(LN2PN_slow_csvName)
		check(err)
		defer LN2PN_slow_csv.Close()
		reader := csv.NewReader(LN2PN_slow_csv)
		fmt.Println("LN2PN_slow_mat")
		for i = 0; i < LN_number; i++ {
			record, err := reader.Read()
			for j = 0; float64(j) < math.Min(float64(PN_number), float64(len(record))) && err == nil; j++ {
				LN2PN_slow_mat[i][j], err = strconv.ParseFloat(record[j], 64)
				check(err)
			}
		}
	} else {
		fmt.Println(LN2PN_slow_csvName, " NOT found!!!")
	}
	// ...
	if exist_file(LN2PN_GABA_csvName) {
		LN2PN_GABA_csv, err = os.Open(LN2PN_GABA_csvName)
		check(err)
		defer LN2PN_GABA_csv.Close()
		reader := csv.NewReader(LN2PN_GABA_csv)
		fmt.Println("LN2PN_GABA_mat")
		for i = 0; i < LN_number; i++ {
			record, err := reader.Read()
			for j = 0; float64(j) < math.Min(float64(PN_number), float64(len(record))) && err == nil; j++ {
				LN2PN_GABA_mat[i][j], err = strconv.ParseFloat(record[j], 64)
				check(err)
			}
		}
	} else {
		fmt.Println(LN2PN_GABA_csvName, " NOT found!!!")
	}
	// ...
	if exist_file(LN2LN_GABA_csvName) {
		LN2LN_GABA_csv, err = os.Open(LN2LN_GABA_csvName)
		check(err)
		defer LN2LN_GABA_csv.Close()
		reader := csv.NewReader(LN2LN_GABA_csv)
		fmt.Println("LN2LN_GABA_mat")
		for i = 0; i < LN_number; i++ {
			record, err := reader.Read()
			for j = 0; float64(j) < math.Min(float64(LN_number), float64(len(record))) && err == nil; j++ {
				LN2LN_GABA_mat[i][j], err = strconv.ParseFloat(record[j], 64)
				check(err)
			}
		}
	} else {
		fmt.Println(LN2LN_GABA_csvName, " NOT found!!!")
	}
	// ...
	if exist_file(PN2PN_nACH_csvName) {
		PN2PN_nACH_csv, err = os.Open(PN2PN_nACH_csvName)
		check(err)
		defer PN2PN_nACH_csv.Close()
		reader := csv.NewReader(PN2PN_nACH_csv)
		fmt.Println("PN2PN_nACH_mat")
		for i = 0; i < PN_number; i++ {
			record, err := reader.Read()
			for j = 0; float64(j) < math.Min(float64(PN_number), float64(len(record))) && err == nil; j++ {
				PN2PN_nACH_mat[i][j], err = strconv.ParseFloat(record[j], 64)
				check(err)
			}
		}
	} else {
		fmt.Println(PN2PN_nACH_csvName, " NOT found!!!")
	}
	// ...
	if exist_file(PN2LN_nACH_csvName) {
		PN2LN_nACH_csv, err = os.Open(PN2LN_nACH_csvName)
		check(err)
		defer PN2LN_nACH_csv.Close()
		reader := csv.NewReader(PN2LN_nACH_csv)
		fmt.Println("PN2LN_nACH_mat")
		for i = 0; i < PN_number; i++ {
			record, err := reader.Read()
			for j = 0; float64(j) < math.Min(float64(LN_number), float64(len(record))) && err == nil; j++ {
				PN2LN_nACH_mat[i][j], err = strconv.ParseFloat(record[j], 64)
				check(err)
			}
		}
	} else {
		fmt.Println(PN2LN_nACH_csvName, " NOT found!!!")
	}
}

func init_matrix(yaml_obj map[interface{}]interface{}) {
	var i, j int64
	defer Exit(Enter("$FN(%v)", yaml_obj))
	if if_readin_matrix > 0 {
		fmt.Println("Using the readin coupling matrix!!!")
		readin_matrix()
		return
	} else {
		fmt.Println("Using randomly set coupling matrix!!!")
	}
	if_slowGABA_overlap = int64(yaml_obj["if_slowGABA_overlap"].(int))
	LN2PN_slow_prob = yaml_obj["LN2PN_slow_prob"].(float64)
	LN2PN_GABA_prob = yaml_obj["LN2PN_GABA_prob"].(float64)
	LN2LN_GABA_prob = yaml_obj["LN2LN_GABA_prob"].(float64)
	PN2PN_nACH_prob = yaml_obj["PN2PN_nACH_prob"].(float64)
	PN2LN_nACH_prob = yaml_obj["PN2LN_nACH_prob"].(float64)
	// ...
	rand.Seed(matrix_seed)
	for i = 0; i < LN_number; i++ { // ln to pn ; slow
		for j = 0; j < PN_number; j++ {
			if rand.Float64() < LN2PN_slow_prob && i != j {
				LN2PN_slow_mat[i][j] = 1
			}
		}
	}
	// ...
	if if_slowGABA_overlap > 0 {
		rand.Seed(matrix_seed) // reseed to ensure LN2PN_slow == LN2PN_GABA
	}
	for i = 0; i < LN_number; i++ { // ln to pn ; GABA
		for j = 0; j < PN_number; j++ {
			if rand.Float64() < LN2PN_GABA_prob && i != j {
				LN2PN_GABA_mat[i][j] = 1
			}
		}
	}
	// ...
	for i = 0; i < LN_number; i++ { // ln to ln
		for j = 0; j < LN_number; j++ {
			if rand.Float64() < LN2LN_GABA_prob && i != j {
				LN2LN_GABA_mat[i][j] = 1
			}
		}
	}
	// ...
	for i = 0; i < PN_number; i++ { // pn to pn
		for j = 0; j < PN_number; j++ {
			if rand.Float64() < PN2PN_nACH_prob && i != j {
				PN2PN_nACH_mat[i][j] = 1
			}
		}
	}
	// ...
	for i = 0; i < PN_number; i++ { // pn to ln
		for j = 0; j < LN_number; j++ {
			if rand.Float64() < PN2LN_nACH_prob && i != j {
				PN2LN_nACH_mat[i][j] = 1
			}
		}
	}
}

func disp_matrix() {
	var i, j int64
	defer Exit(Enter("$FN()"))
	fmt.Println("\nLN2PN_slow_mat")
	for i = 0; i < LN_number; i++ { // ln to pn ; slow
		for j = 0; j < PN_number; j++ {
			fmt.Print(LN2PN_slow_mat[i][j], " ")
		}
		fmt.Println()
	}
	// ...
	fmt.Println("\nLN2PN_GABA_mat")
	for i = 0; i < LN_number; i++ { // ln to pn ; GABA
		for j = 0; j < PN_number; j++ {
			fmt.Print(LN2PN_GABA_mat[i][j], " ")
		}
		fmt.Println()
	}
	// ...
	fmt.Println("\nLN2LN_GABA_mat")
	for i = 0; i < LN_number; i++ { // ln to ln
		for j = 0; j < LN_number; j++ {
			fmt.Print(LN2LN_GABA_mat[i][j], " ")
		}
		fmt.Println()
	}
	// ...
	fmt.Println("\nPN2PN_nACH_mat")
	for i = 0; i < PN_number; i++ { // pn to pn
		for j = 0; j < PN_number; j++ {
			fmt.Print(PN2PN_nACH_mat[i][j], " ")
		}
		fmt.Println()
	}
	// ...
	fmt.Println("\nPN2LN_nACH_mat")
	for i = 0; i < PN_number; i++ { // pn to ln
		for j = 0; j < LN_number; j++ {
			fmt.Print(PN2LN_nACH_mat[i][j], " ")
		}
		fmt.Println()
	}
	fmt.Println()
}

func save_matrix() {
	var i, j int64
	defer Exit(Enter("$FN()"))
	LN2PN_slow_txt, err = os.OpenFile(LN2PN_slow_txtName, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	LN2PN_GABA_txt, err = os.OpenFile(LN2PN_GABA_txtName, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	LN2LN_GABA_txt, err = os.OpenFile(LN2LN_GABA_txtName, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	PN2PN_nACH_txt, err = os.OpenFile(PN2PN_nACH_txtName, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	PN2LN_nACH_txt, err = os.OpenFile(PN2LN_nACH_txtName, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	// ...
	for i = 0; i < LN_number; i++ { // ln to pn ; slow
		for j = 0; j < PN_number; j++ {
			_, _ = LN2PN_slow_txt.Write([]byte(strconv.FormatFloat(LN2PN_slow_mat[i][j], 'f', 3, 64)))
			_, _ = LN2PN_slow_txt.WriteString(" ")
		}
		_, _ = LN2PN_slow_txt.WriteString("\n")
	}
	// ...
	for i = 0; i < LN_number; i++ { // ln to pn ; GABA
		for j = 0; j < PN_number; j++ {
			_, _ = LN2PN_GABA_txt.Write([]byte(strconv.FormatFloat(LN2PN_GABA_mat[i][j], 'f', 3, 64)))
			_, _ = LN2PN_GABA_txt.WriteString(" ")
		}
		_, _ = LN2PN_GABA_txt.WriteString("\n")
	}
	// ...
	for i = 0; i < LN_number; i++ { // ln to ln
		for j = 0; j < LN_number; j++ {
			_, _ = LN2LN_GABA_txt.Write([]byte(strconv.FormatFloat(LN2LN_GABA_mat[i][j], 'f', 3, 64)))
			_, _ = LN2LN_GABA_txt.WriteString(" ")
		}
		_, _ = LN2LN_GABA_txt.WriteString("\n")
	}
	// ...
	for i = 0; i < PN_number; i++ { // pn to pn
		for j = 0; j < PN_number; j++ {
			_, _ = PN2PN_nACH_txt.Write([]byte(strconv.FormatFloat(PN2PN_nACH_mat[i][j], 'f', 3, 64)))
			_, _ = PN2PN_nACH_txt.WriteString(" ")
		}
		_, _ = PN2PN_nACH_txt.WriteString("\n")
	}
	// ...
	for i = 0; i < PN_number; i++ { // pn to ln
		for j = 0; j < LN_number; j++ {
			_, _ = PN2LN_nACH_txt.Write([]byte(strconv.FormatFloat(PN2LN_nACH_mat[i][j], 'f', 3, 64)))
			_, _ = PN2LN_nACH_txt.WriteString(" ")
		}
		_, _ = PN2LN_nACH_txt.WriteString("\n")
	}
	// ...
	LN2PN_slow_txt.Close()
	LN2PN_GABA_txt.Close()
	LN2LN_GABA_txt.Close()
	PN2PN_nACH_txt.Close()
	PN2LN_nACH_txt.Close()
	fmt.Println("Coupling matrixes saved to files.")
	fmt.Println()
}

func forward_states_PN(id int64) {
	var clock float64 = float64(int64(math.Floor(CLOCK)) % period_len)
	defer Exit(Enter("$FN(%v)", id))
	if id < 0 || id >= PN_number {
		return
	}
	PN_reactor[id][1].m_K = PN_reactor[id][0].m_K + PN_reactor[id][0].dfdt_m_K()*time_stepLen
	PN_reactor[id][1].h_K = PN_reactor[id][0].h_K + PN_reactor[id][0].dfdt_h_K()*time_stepLen
	PN_reactor[id][1].m_Na = PN_reactor[id][0].m_Na + PN_reactor[id][0].dfdt_m_Na()*time_stepLen
	PN_reactor[id][1].h_Na = PN_reactor[id][0].h_Na + PN_reactor[id][0].dfdt_h_Na()*time_stepLen
	PN_reactor[id][1].m_A = PN_reactor[id][0].m_A + PN_reactor[id][0].dfdt_m_A()*time_stepLen
	PN_reactor[id][1].h_A = PN_reactor[id][0].h_A + PN_reactor[id][0].dfdt_h_A()*time_stepLen
	PN_reactor[id][1].O_nACH = PN_reactor[id][0].O_nACH + PN_reactor[id][0].dfdt_O_nACH()*time_stepLen
	PN_reactor[id][1].V = PN_reactor[id][0].V + PN_reactor[id][0].dfdt_V()*time_stepLen
	// ...
	if PN_reactor[id][0].V <= V_threshold_PN && PN_reactor[id][1].V > V_threshold_PN {
		PN_reactor[id][1].lastSpike = CLOCK
		if freq_scoping == "noStim" && (clock < stim_onset || clock >= (stim_offset+stim_decay)) {
			PN_spike_freq[id] += float64(ms_per_second) / float64(stim_onset+float64(period_len)-stim_offset-stim_decay)
		} else if freq_scoping == "rising" && clock >= stim_onset && clock < (stim_onset+stim_rise) {
			PN_spike_freq[id] += float64(ms_per_second) / float64(stim_rise) // plus 1 in the form freq in S
		} else if freq_scoping == "stimulated" && clock >= (stim_onset+stim_rise) && clock < stim_offset {
			PN_spike_freq[id] += float64(ms_per_second) / float64(stim_offset-stim_onset-stim_rise) // plus 1 in the form freq in S
		} else if freq_scoping == "decay" && clock >= stim_offset && clock < (stim_offset+stim_decay) {
			PN_spike_freq[id] += float64(ms_per_second) / float64(stim_decay) // plus 1 in the form freq in S
		}
	} else {
		PN_reactor[id][1].lastSpike = PN_reactor[id][0].lastSpike
	}
}

func forward_states_LN(id int64) {
	var clock float64 = float64(int64(math.Floor(CLOCK)) % period_len)
	defer Exit(Enter("$FN(%v)", id))
	if id < 0 || id >= LN_number {
		return
	}
	LN_reactor[id][1].m_K = LN_reactor[id][0].m_K + LN_reactor[id][0].dfdt_m_K()*time_stepLen
	LN_reactor[id][1].h_K = LN_reactor[id][0].h_K + LN_reactor[id][0].dfdt_h_K()*time_stepLen
	LN_reactor[id][1].m_Ca = LN_reactor[id][0].m_Ca + LN_reactor[id][0].dfdt_m_Ca()*time_stepLen
	LN_reactor[id][1].h_Ca = LN_reactor[id][0].h_Ca + LN_reactor[id][0].dfdt_h_Ca()*time_stepLen
	LN_reactor[id][1].m_B = LN_reactor[id][0].m_B + LN_reactor[id][0].dfdt_m_B()*time_stepLen
	LN_reactor[id][1].h_B = LN_reactor[id][0].h_B + LN_reactor[id][0].dfdt_h_B()*time_stepLen
	LN_reactor[id][1].F_Ca = LN_reactor[id][0].F_Ca + LN_reactor[id][0].dfdt_F_Ca()*time_stepLen
	LN_reactor[id][1].O_GABA = LN_reactor[id][0].O_GABA + LN_reactor[id][0].dfdt_O_GABA()*time_stepLen
	LN_reactor[id][1].G_slow = LN_reactor[id][0].G_slow + LN_reactor[id][0].dfdt_G_slow()*time_stepLen
	LN_reactor[id][1].R_slow = LN_reactor[id][0].R_slow + LN_reactor[id][0].dfdt_R_slow()*time_stepLen
	LN_reactor[id][1].V = LN_reactor[id][0].V + LN_reactor[id][0].dfdt_V()*time_stepLen // mereged after [new2]
	// ...
	if LN_reactor[id][0].V <= V_threshold_LN && LN_reactor[id][1].V > V_threshold_LN { // is this set properly?
		LN_reactor[id][1].lastSpike = CLOCK
		if freq_scoping == "noStim" && (clock < stim_onset || clock >= (stim_offset+stim_decay)) {
			LN_spike_freq[id] += float64(ms_per_second) / float64(stim_onset+float64(period_len)-stim_offset-stim_decay)
		} else if freq_scoping == "rising" && clock >= stim_onset && clock < (stim_onset+stim_rise) {
			LN_spike_freq[id] += float64(ms_per_second) / float64(stim_rise) // plus 1 in the form freq in S
		} else if freq_scoping == "stimulated" && clock >= (stim_onset+stim_rise) && clock < stim_offset {
			LN_spike_freq[id] += float64(ms_per_second) / float64(stim_offset-stim_onset-stim_rise) // plus 1 in the form freq in S
		} else if freq_scoping == "decay" && clock >= stim_offset && clock < (stim_offset+stim_decay) {
			LN_spike_freq[id] += float64(ms_per_second) / float64(stim_decay) // plus 1 in the form freq in S
		}
	} else { // reveised after [try]
		LN_reactor[id][1].lastSpike = LN_reactor[id][0].lastSpike
	}
}

func radiate_currents_PN(id int64) {
	var ToGABA, TonACH, TgSlow float64 = 0, 0, 0
	var i int64
	defer Exit(Enter("$FN(%v)", id))
	if id < 0 || id >= PN_number {
		return
	}
	for i = 0; i < LN_number; i++ { // ln 2 pn
		ToGABA += LN2PN_GABA_mat[i][id] * LN_reactor[i][0].O_GABA
		TgSlow += LN2PN_slow_mat[i][id] * LN_reactor[i][0].G_slow
	} // yes, there should times xN_reactor[i][0].x; kind of update states
	for i = 0; i < PN_number; i++ { // pn 2 pn
		TonACH += PN2PN_nACH_mat[i][id] * PN_reactor[i][0].O_nACH // how to process the matrix here
	}
	// ...
	PN_reactor[id][1].I_K = g_K_PN * math.Pow(PN_reactor[id][1].m_K, 4) * math.Pow(PN_reactor[id][1].h_K, 0) * (PN_reactor[id][1].V + E_K_PN)
	PN_reactor[id][1].I_Na = g_Na_PN * math.Pow(PN_reactor[id][1].m_Na, 3) * math.Pow(PN_reactor[id][1].h_Na, 1) * (PN_reactor[id][1].V + E_Na_PN)
	PN_reactor[id][1].I_A = g_A_PN * math.Pow(PN_reactor[id][1].m_A, 4) * math.Pow(PN_reactor[id][1].h_A, 1) * (PN_reactor[id][1].V + E_A_PN)
	PN_reactor[id][1].I_GABA = g_GABA_PN * ToGABA * (PN_reactor[id][1].V + E_GABA_PN)
	PN_reactor[id][1].I_nACH = g_nACH_PN * TonACH * (PN_reactor[id][1].V + E_nACH_PN)
	PN_reactor[id][1].I_slow = g_slow_PN * (math.Pow(TgSlow, 4) / (math.Pow(TgSlow, 4) + 100.0)) * (PN_reactor[id][1].V + E_slow_PN)
}

func radiate_currents_LN(id int64) {
	var ToGABA, TonACH float64 = 0, 0
	var i int64
	defer Exit(Enter("$FN(%v)", id))
	if id < 0 || id >= LN_number {
		return
	}
	for i = 0; i < LN_number; i++ { // how to process this thing LN2LN
		ToGABA += LN2LN_GABA_mat[i][id] * LN_reactor[i][0].O_GABA
	} // yes, there should times xN_reactor[i][0].x; kind of update states
	for i = 0; i < PN_number; i++ { // PN 2 LN
		TonACH += PN2LN_nACH_mat[i][id] * PN_reactor[i][0].O_nACH
	}
	// ...
	LN_reactor[id][1].I_K = g_K_LN * math.Pow(LN_reactor[id][1].m_K, 4) * math.Pow(LN_reactor[id][1].h_K, 0) * (LN_reactor[id][1].V + E_K_LN)
	LN_reactor[id][1].I_Ca = g_Ca_LN * math.Pow(LN_reactor[id][1].m_Ca, 2) * math.Pow(LN_reactor[id][1].h_Ca, 1) * (LN_reactor[id][1].V + E_Ca_LN)
	LN_reactor[id][1].I_B = g_B_LN * math.Pow(LN_reactor[id][1].m_B, 1) * math.Pow(LN_reactor[id][1].h_B, 0) * (LN_reactor[id][1].V + E_B_LN)
	LN_reactor[id][1].I_GABA = g_GABA_LN * ToGABA * (LN_reactor[id][1].V + E_GABA_LN)
	LN_reactor[id][1].I_nACH = g_nACH_LN * TonACH * (LN_reactor[id][1].V + E_nACH_LN)
}

func get_stim_PN(id int64) {
	var inputRate, ret float64 = 0, fading_rate * PN_reactor[id][0].stim
	var clock float64 = float64(int64(math.Floor(CLOCK)) % period_len)
	var i int64 = 0
	defer Exit(Enter("$FN(%v)", id))
	if id < 0 || id >= PN_number {
		return
	}
	// ...
	if rand.Float64() < (BG_input_rate * time_stepLen) {
		ret -= BG_input_strength_PN
	}
	if clock < stim_onset {
		inputRate = 0
	} else if clock >= stim_onset && clock < (stim_onset+stim_rise) { // rising stimulus
		//inputRate = math.Exp(-math.Pow(clock-stim_onset-stim_rise, 2)/1000000) * ORN_input_rate
		inputRate = ((clock - stim_onset) / stim_rise) * ORN_input_rate
	} else if clock >= (stim_onset+stim_rise) && clock < stim_offset { // normal stimulus
		inputRate = ORN_input_rate
	} else if clock >= stim_offset && clock < (stim_offset+stim_decay) { // decay stimulus
		//inputRate = math.Exp(-math.Sqrt(clock-stim_offset)/math.Sqrt(1000)) * ORN_input_rate
		inputRate = ((stim_offset + stim_decay - clock) / stim_decay) * ORN_input_rate
	} else if clock >= (stim_offset + stim_decay) { // after decay
		inputRate = 0
	}
	if odor_PN_set.Contains(id) {
		for i = 0; i < ORN_number; i++ {
			if rand.Float64() < (inputRate * time_stepLen) {
				ret -= ORN_input_strength_PN
			}
		}
	}
	PN_reactor[id][1].stim = ret
	// ...
	if if_random_stim > 0 {
		return // Poisson stimulus or steady stimulus
	} else if clock < stim_onset {
		PN_reactor[id][1].stim = PN_stim_baseline
	} else if clock >= stim_onset && clock < (stim_onset+stim_rise) { // rising stimulus
		PN_reactor[id][1].stim = PN_stim_baseline + ((clock-stim_onset)/stim_rise)*(PN_stim_plusORNs-PN_stim_baseline)
	} else if clock >= (stim_onset+stim_rise) && clock < stim_offset { // normal stimulus
		PN_reactor[id][1].stim = PN_stim_plusORNs
	} else if clock >= stim_offset && clock < (stim_offset+stim_decay) { // decay stimulus
		PN_reactor[id][1].stim = PN_stim_baseline + ((stim_offset+stim_decay-clock)/stim_decay)*(PN_stim_plusORNs-PN_stim_baseline)
	} else if clock >= (stim_offset + stim_decay) { // after decay
		PN_reactor[id][1].stim = PN_stim_baseline
	}
}

func get_stim_LN(id int64) {
	var inputRate, ret float64 = 0, fading_rate * LN_reactor[id][0].stim
	var clock float64 = float64(int64(math.Floor(CLOCK)) % period_len)
	var i int64 = 0
	defer Exit(Enter("$FN(%v)", id))
	if id < 0 || id >= LN_number {
		return
	}
	// ...
	if rand.Float64() < (BG_input_rate * time_stepLen) {
		ret -= BG_input_strength_LN
	}
	if clock < stim_onset {
		inputRate = 0
	} else if clock >= stim_onset && clock < (stim_onset+stim_rise) { // rising stimulus
		//inputRate = math.Exp(-math.Pow(clock-stim_onset-stim_rise, 2)/1000000) * ORN_input_rate
		inputRate = ((clock - stim_onset) / stim_rise) * ORN_input_rate
	} else if clock >= (stim_onset+stim_rise) && clock < stim_offset { // normal stimulus
		inputRate = ORN_input_rate
	} else if clock >= stim_offset && clock < (stim_offset+stim_decay) { // decay stimulus
		//inputRate = math.Exp(-math.Sqrt(clock-stim_offset)/math.Sqrt(1000)) * ORN_input_rate
		inputRate = ((stim_offset + stim_decay - clock) / stim_decay) * ORN_input_rate
	} else if clock >= (stim_offset + stim_decay) { // after decay
		inputRate = 0
	}
	if odor_LN_set.Contains(id) {
		for i = 1; i <= ORN_number; i++ {
			if rand.Float64() < (inputRate * time_stepLen) {
				ret -= ORN_input_strength_LN
			}
		}
	}
	LN_reactor[id][1].stim = ret
	// ...
	if if_random_stim > 0 {
		return // Poisson stimulus or steady stimulus
	} else if clock < stim_onset {
		LN_reactor[id][1].stim = LN_stim_baseline
	} else if clock >= stim_onset && clock < (stim_onset+stim_rise) { // rising stimulus
		LN_reactor[id][1].stim = LN_stim_baseline + ((clock-stim_onset)/stim_rise)*(LN_stim_plusORNs-LN_stim_baseline)
	} else if clock >= (stim_onset+stim_rise) && clock < stim_offset { // normal stimulus
		LN_reactor[id][1].stim = LN_stim_plusORNs
	} else if clock >= stim_offset && clock < (stim_offset+stim_decay) { // decay stimulus
		LN_reactor[id][1].stim = LN_stim_baseline + ((stim_offset+stim_decay-clock)/stim_decay)*(LN_stim_plusORNs-LN_stim_baseline)
	} else if clock >= (stim_offset + stim_decay) { // after decay
		LN_reactor[id][1].stim = LN_stim_baseline
	}
}

func shift_reactor_PN(id, click int64) {
	defer Exit(Enter("$FN(%v, %v)", id, click))
	if id < 0 || id >= PN_number {
		return
	}
	if click%100 == 0 { // 1ms 1 record
		PN_stim_list[click/100][id] = PN_reactor[id][1].stim // at this click, we input this stim, and
		PN_V_list[click/100][id] = PN_reactor[id][1].V       // at next click, we will get the corresponding V
		PN_slow_list[click/100][id] = PN_reactor[id][1].I_slow
		PN_GABA_list[click/100][id] = PN_reactor[id][1].I_GABA
		PN_nACH_list[click/100][id] = PN_reactor[id][1].I_nACH
	}
	// ...
	PN_reactor[id][0] = PN_reactor[id][1]
	PN_reactor[id][1].m_K = 0
	PN_reactor[id][1].h_K = 1
	PN_reactor[id][1].m_Na = 0
	PN_reactor[id][1].h_Na = 1
	PN_reactor[id][1].m_A = 0
	PN_reactor[id][1].h_A = 1
	PN_reactor[id][1].O_nACH = 0
	PN_reactor[id][1].V = 0
	PN_reactor[id][1].I_K = 0
	PN_reactor[id][1].I_Na = 0
	PN_reactor[id][1].I_A = 0
	PN_reactor[id][1].I_GABA = 0
	PN_reactor[id][1].I_nACH = 0
	PN_reactor[id][1].I_slow = 0
}

func shift_reactor_LN(id, click int64) {
	defer Exit(Enter("$FN(%v, %v)", id, click))
	if id < 0 || id >= LN_number {
		return
	}
	if click%100 == 0 { // 1ms 1 record
		LN_stim_list[click/100][id] = LN_reactor[id][1].stim // at this click, we input this stim, and
		LN_V_list[click/100][id] = LN_reactor[id][1].V       // at next click, we will get the corresponding V
		LN_GABA_list[click/100][id] = LN_reactor[id][1].I_GABA
		LN_nACH_list[click/100][id] = LN_reactor[id][1].I_nACH
	}
	// ...
	LN_reactor[id][0] = LN_reactor[id][1]
	LN_reactor[id][1].m_K = 0
	LN_reactor[id][1].h_K = 1
	LN_reactor[id][1].m_Ca = 0
	LN_reactor[id][1].h_Ca = 1
	LN_reactor[id][1].m_B = 0
	LN_reactor[id][1].h_B = 1
	LN_reactor[id][1].F_Ca = 0
	LN_reactor[id][1].O_GABA = 0
	LN_reactor[id][1].G_slow = 0
	LN_reactor[id][1].R_slow = 0
	LN_reactor[id][1].V = 0
	LN_reactor[id][1].I_K = 0
	LN_reactor[id][1].I_Ca = 0
	LN_reactor[id][1].I_B = 0
	LN_reactor[id][1].I_GABA = 0
	LN_reactor[id][1].I_nACH = 0
}

func take_iteration(click int64) { // xNs_pre -> xNs_cur
	var id int64
	defer Exit(Enter("$FN(%v)", click))
	for id = 0; id < PN_number; id++ {
		forward_states_PN(id)
		radiate_currents_PN(id)
		get_stim_PN(id)
		shift_reactor_PN(id, click)
	}
	// ...
	for id = 0; id < LN_number; id++ {
		forward_states_LN(id)
		radiate_currents_LN(id)
		get_stim_LN(id)
		shift_reactor_LN(id, click)
	}
}

////////////////////////////////////////////////////////////
// For this part``` ////////////////////////////////////////
///////////////// do the input, interface and main func ////
////////////////////////////////////////////////////////////

func check(e error) {
	if e != nil {
		panic(e)
	}
}

func exist_file(filename string) bool {
	_, err := os.Stat(filename)
	return err == nil || os.IsExist(err)
}

var (
	ms_per_frame int64 = 750 // millisecond per gnuplot frame
	// ...
	filename_config             string = "config.yaml"
	filename_vRaster_PN         string = strings.Join([]string{plot_file_dir, "plot_vRaster_PN.txt"}, "")    // will be added an appendix
	filename_vRaster_LN         string = strings.Join([]string{plot_file_dir, "plot_vRaster_LN.txt"}, "")    // in proc_filenames(),
	filename_stimRaster_PN      string = strings.Join([]string{plot_file_dir, "plot_stimRaster_PN.txt"}, "") // so that
	filename_stimRaster_LN      string = strings.Join([]string{plot_file_dir, "plot_stimRaster_LN.txt"}, "") // many processes
	filename_synaRaster_slow_PN string = strings.Join([]string{plot_file_dir, "plot_synaRaster_slow_PN.txt"}, "")
	filename_synaRaster_GABA_PN string = strings.Join([]string{plot_file_dir, "plot_synaRaster_GABA_PN.txt"}, "")
	filename_synaRaster_GABA_LN string = strings.Join([]string{plot_file_dir, "plot_synaRaster_GABA_LN.txt"}, "")
	filename_synaRaster_nACH_PN string = strings.Join([]string{plot_file_dir, "plot_synaRaster_nACH_PN.txt"}, "")
	filename_synaRaster_nACH_LN string = strings.Join([]string{plot_file_dir, "plot_synaRaster_nACH_LN.txt"}, "")
	// ...
	filename_vFluct_PN         string = strings.Join([]string{plot_file_dir, "plot_vFluct_PN.txt"}, "") // can be run
	filename_vFluct_LN         string = strings.Join([]string{plot_file_dir, "plot_vFluct_LN.txt"}, "") // simultaneously...
	filename_stimFluct_PN      string = strings.Join([]string{plot_file_dir, "plot_stimFluct_PN.txt"}, "")
	filename_stimFluct_LN      string = strings.Join([]string{plot_file_dir, "plot_stimFluct_LN.txt"}, "")
	filename_synaFluct_slow_PN string = strings.Join([]string{plot_file_dir, "plot_synaFluct_slow_PN.txt"}, "")
	filename_synaFluct_GABA_PN string = strings.Join([]string{plot_file_dir, "plot_synaFluct_GABA_PN.txt"}, "")
	filename_synaFluct_GABA_LN string = strings.Join([]string{plot_file_dir, "plot_synaFluct_GABA_LN.txt"}, "")
	filename_synaFluct_nACH_PN string = strings.Join([]string{plot_file_dir, "plot_synaFluct_nACH_PN.txt"}, "")
	filename_synaFluct_nACH_LN string = strings.Join([]string{plot_file_dir, "plot_synaFluct_nACH_LN.txt"}, "")
	// ...
	file_config                            *os.File
	file_vRaster_PN, file_vRaster_LN       *os.File
	plot_vRaster_PN, plot_vRaster_LN       *gnuplot.Plotter
	file_stimRaster_PN, file_stimRaster_LN *os.File
	plot_stimRaster_PN, plot_stimRaster_LN *gnuplot.Plotter
	file_vFluct_PN, file_vFluct_LN         *os.File
	plot_vFluct_PN, plot_vFluct_LN         *gnuplot.Plotter
	file_stimFluct_PN, file_stimFluct_LN   *os.File
	plot_stimFluct_PN, plot_stimFluct_LN   *gnuplot.Plotter
	// ...
	file_synaRaster_slow_PN                          *os.File
	plot_synaRaster_slow_PN                          *gnuplot.Plotter
	file_synaRaster_GABA_PN, file_synaRaster_GABA_LN *os.File
	plot_synaRaster_GABA_PN, plot_synaRaster_GABA_LN *gnuplot.Plotter
	file_synaRaster_nACH_PN, file_synaRaster_nACH_LN *os.File
	plot_synaRaster_nACH_PN, plot_synaRaster_nACH_LN *gnuplot.Plotter
	// ...
	file_synaFluct_slow_PN                         *os.File
	plot_synaFluct_slow_PN                         *gnuplot.Plotter
	file_synaFluct_GABA_PN, file_synaFluct_GABA_LN *os.File
	plot_synaFluct_GABA_PN, plot_synaFluct_GABA_LN *gnuplot.Plotter
	file_synaFluct_nACH_PN, file_synaFluct_nACH_LN *os.File
	plot_synaFluct_nACH_PN, plot_synaFluct_nACH_LN *gnuplot.Plotter
	// save data for further processing...
	doc_V_PN, doc_V_LN       *os.File
	doc_stim_PN, doc_stim_LN *os.File
	docname_V_PN             string = strings.Join([]string{save_file_dir, "doc_V_PN.txt"}, "")
	docname_V_LN             string = strings.Join([]string{save_file_dir, "doc_V_LN.txt"}, "")
	docname_stim_PN          string = strings.Join([]string{save_file_dir, "doc_stim_PN.txt"}, "")
	docname_stim_LN          string = strings.Join([]string{save_file_dir, "doc_stim_LN.txt"}, "")
	// ...
	doc_slow_PN              *os.File
	doc_GABA_PN, doc_GABA_LN *os.File
	doc_nACH_PN, doc_nACH_LN *os.File
	docname_slow_PN          string = strings.Join([]string{save_file_dir, "doc_stim_PN.txt"}, "")
	docname_GABA_PN          string = strings.Join([]string{save_file_dir, "doc_GABA_PN.txt"}, "")
	docname_GABA_LN          string = strings.Join([]string{save_file_dir, "doc_GABA_LN.txt"}, "")
	docname_nACH_PN          string = strings.Join([]string{save_file_dir, "doc_nACH_PN.txt"}, "")
	docname_nACH_LN          string = strings.Join([]string{save_file_dir, "doc_nACH_LN.txt"}, "")
	// ...
	LN2PN_slow_txt     *os.File
	LN2PN_GABA_txt     *os.File
	LN2LN_GABA_txt     *os.File
	PN2PN_nACH_txt     *os.File
	PN2LN_nACH_txt     *os.File
	LN2PN_slow_txtName string = strings.Join([]string{save_file_dir, "mat_LN2PN_slow.txt"}, "")
	LN2PN_GABA_txtName string = strings.Join([]string{save_file_dir, "mat_LN2PN_GABA.txt"}, "")
	LN2LN_GABA_txtName string = strings.Join([]string{save_file_dir, "mat_LN2LN_GABA.txt"}, "")
	PN2PN_nACH_txtName string = strings.Join([]string{save_file_dir, "mat_PN2PN_nACH.txt"}, "")
	PN2LN_nACH_txtName string = strings.Join([]string{save_file_dir, "mat_PN2LN_nACH.txt"}, "")
	// ...
	PN_spike_freq_file     *os.File
	LN_spike_freq_file     *os.File
	PN_spike_freq_filename string = strings.Join([]string{save_file_dir, "PN_spike_freq.txt"}, "")
	LN_spike_freq_filename string = strings.Join([]string{save_file_dir, "LN_spike_freq.txt"}, "")
	// ...
	LN2PN_slow_csv     *os.File
	LN2PN_GABA_csv     *os.File
	LN2LN_GABA_csv     *os.File
	PN2PN_nACH_csv     *os.File
	PN2LN_nACH_csv     *os.File
	LN2PN_slow_csvName string = strings.Join([]string{save_file_dir, "mat_LN2PN_slow.csv"}, "")
	LN2PN_GABA_csvName string = strings.Join([]string{save_file_dir, "mat_LN2PN_GABA.csv"}, "")
	LN2LN_GABA_csvName string = strings.Join([]string{save_file_dir, "mat_LN2LN_GABA.csv"}, "")
	PN2PN_nACH_csvName string = strings.Join([]string{save_file_dir, "mat_PN2PN_nACH.csv"}, "")
	PN2LN_nACH_csvName string = strings.Join([]string{save_file_dir, "mat_PN2LN_nACH.csv"}, "")
)

func init_rand() {
	defer Exit(Enter("$FN()"))
	rand.Seed(time.Now().Unix())
}

func init_para() { // add more thing to write at initial?
	defer Exit(Enter("$FN()"))
	if exist_file(filename_config) && use_exist_config {
		fmt.Println("WARNING: using old config file!!!")
	} else {
		fmt.Println("WARNING: creat new config file!!!")
		file_config, err = os.OpenFile(filename_config, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
		check(err) // always rewrite the config.yaml !!! !!! !!!
		// ...
		_, _ = file_config.WriteString("if_running: ")
		_, _ = file_config.Write([]byte(strconv.FormatInt(if_running, 10)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("if_slowGABA_overlap: ")
		_, _ = file_config.Write([]byte(strconv.FormatInt(if_slowGABA_overlap, 10)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("if_random_stim: ")
		_, _ = file_config.Write([]byte(strconv.FormatInt(if_random_stim, 10)))
		_, _ = file_config.WriteString("\n\n")
		// ...
		_, _ = file_config.WriteString("if_realtime_plot: ")
		_, _ = file_config.Write([]byte(strconv.FormatInt(if_realtime_plot, 10)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("if_PN_plot: ")
		_, _ = file_config.Write([]byte(strconv.FormatInt(if_PN_plot, 10)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("if_LN_plot: ")
		_, _ = file_config.Write([]byte(strconv.FormatInt(if_LN_plot, 10)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("if_stimRaster_plot: ")
		_, _ = file_config.Write([]byte(strconv.FormatInt(if_stimRaster_plot, 10)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("if_vRaster_plot: ")
		_, _ = file_config.Write([]byte(strconv.FormatInt(if_vRaster_plot, 10)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("if_synaRaster_plot: ")
		_, _ = file_config.Write([]byte(strconv.FormatInt(if_synaRaster_plot, 10)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("if_vFluct_plot: ")
		_, _ = file_config.Write([]byte(strconv.FormatInt(if_vFluct_plot, 10)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("if_stimFluct_plot: ")
		_, _ = file_config.Write([]byte(strconv.FormatInt(if_stimFluct_plot, 10)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("if_synaFluct_plot: ")
		_, _ = file_config.Write([]byte(strconv.FormatInt(if_synaFluct_plot, 10)))
		_, _ = file_config.WriteString("\n\n")
		// ...
		_, _ = file_config.WriteString("if_LN2PN_slowed: ")
		_, _ = file_config.Write([]byte(strconv.FormatInt(if_LN2PN_slowed, 10)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("if_LN2PN_GABAed: ")
		_, _ = file_config.Write([]byte(strconv.FormatInt(if_LN2PN_GABAed, 10)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("if_LN2LN_GABAed: ")
		_, _ = file_config.Write([]byte(strconv.FormatInt(if_LN2LN_GABAed, 10)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("if_PN2PN_nACHed: ")
		_, _ = file_config.Write([]byte(strconv.FormatInt(if_PN2PN_nACHed, 10)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("if_PN2LN_nACHed: ")
		_, _ = file_config.Write([]byte(strconv.FormatInt(if_PN2LN_nACHed, 10)))
		_, _ = file_config.WriteString("\n\n")
		// ...
		_, _ = file_config.WriteString("scale_slow_PN: ")
		_, _ = file_config.Write([]byte(strconv.FormatFloat(scale_slow_PN, 'f', 6, 64)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("scale_GABA_PN: ")
		_, _ = file_config.Write([]byte(strconv.FormatFloat(scale_GABA_PN, 'f', 6, 64)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("scale_GABA_LN: ")
		_, _ = file_config.Write([]byte(strconv.FormatFloat(scale_GABA_LN, 'f', 6, 64)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("scale_nACH_PN: ")
		_, _ = file_config.Write([]byte(strconv.FormatFloat(scale_nACH_PN, 'f', 6, 64)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("scale_nACH_LN: ")
		_, _ = file_config.Write([]byte(strconv.FormatFloat(scale_nACH_LN, 'f', 6, 64)))
		_, _ = file_config.WriteString("\n\n")
		// ...
		_, _ = file_config.WriteString("LN2PN_slow_prob: ")
		_, _ = file_config.Write([]byte(strconv.FormatFloat(LN2PN_slow_prob, 'f', 6, 64)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("LN2PN_GABA_prob: ")
		_, _ = file_config.Write([]byte(strconv.FormatFloat(LN2PN_GABA_prob, 'f', 6, 64)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("LN2LN_GABA_prob: ")
		_, _ = file_config.Write([]byte(strconv.FormatFloat(LN2LN_GABA_prob, 'f', 6, 64)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("PN2PN_nACH_prob: ")
		_, _ = file_config.Write([]byte(strconv.FormatFloat(PN2PN_nACH_prob, 'f', 6, 64)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("PN2LN_nACH_prob: ")
		_, _ = file_config.Write([]byte(strconv.FormatFloat(PN2LN_nACH_prob, 'f', 6, 64)))
		_, _ = file_config.WriteString("\n\n")
		// ...
		_, _ = file_config.WriteString("BG_input_rate: ")
		_, _ = file_config.Write([]byte(strconv.FormatFloat(BG_input_rate, 'f', 6, 64)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("ORN_input_rate: ")
		_, _ = file_config.Write([]byte(strconv.FormatFloat(ORN_input_rate, 'f', 6, 64)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("BG_input_strength_PN: ")
		_, _ = file_config.Write([]byte(strconv.FormatFloat(BG_input_strength_PN, 'f', 6, 64)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("BG_input_strength_LN: ")
		_, _ = file_config.Write([]byte(strconv.FormatFloat(BG_input_strength_LN, 'f', 6, 64)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("ORN_input_strength_PN: ")
		_, _ = file_config.Write([]byte(strconv.FormatFloat(ORN_input_strength_PN, 'f', 6, 64)))
		_, _ = file_config.WriteString("\n")
		_, _ = file_config.WriteString("ORN_input_strength_LN: ")
		_, _ = file_config.Write([]byte(strconv.FormatFloat(ORN_input_strength_LN, 'f', 6, 64)))
		_, _ = file_config.WriteString("\n\n")
		file_config.Sync()
		file_config.Close()
	}
}

func decode_yaml(mm map[interface{}]interface{}) map[interface{}]interface{} {
	defer Exit(Enter("$FN(%v)", mm))
	data, err := ioutil.ReadFile(filename_config)
	if err != nil { // when yaml file is opened outside, it cannot be opened here and there will be an err
		return mm // change nothing in this case
	} // else if the yaml file opened, then read it
	m := make(map[interface{}]interface{})
	err = yaml.Unmarshal([]byte(data), &m)
	return m
}

func proc_para(yaml_obj map[interface{}]interface{}) {
	defer Exit(Enter("$FN(%v)", yaml_obj))
	if_running = int64(yaml_obj["if_running"].(int))
	if_random_stim = int64(yaml_obj["if_random_stim"].(int))
	if_realtime_plot = int64(yaml_obj["if_realtime_plot"].(int))
	if_PN_plot = int64(yaml_obj["if_PN_plot"].(int))
	if_LN_plot = int64(yaml_obj["if_LN_plot"].(int))
	if_stimRaster_plot = int64(yaml_obj["if_stimRaster_plot"].(int))
	if_vRaster_plot = int64(yaml_obj["if_vRaster_plot"].(int))
	if_synaRaster_plot = int64(yaml_obj["if_synaRaster_plot"].(int))
	if_vFluct_plot = int64(yaml_obj["if_vFluct_plot"].(int))
	if_stimFluct_plot = int64(yaml_obj["if_stimFluct_plot"].(int))
	if_synaFluct_plot = int64(yaml_obj["if_synaFluct_plot"].(int))
	// ...
	if_LN2PN_slowed = int64(yaml_obj["if_LN2PN_slowed"].(int))
	if_LN2PN_GABAed = int64(yaml_obj["if_LN2PN_GABAed"].(int))
	if_LN2LN_GABAed = int64(yaml_obj["if_LN2LN_GABAed"].(int))
	if_PN2PN_nACHed = int64(yaml_obj["if_PN2PN_nACHed"].(int))
	if_PN2LN_nACHed = int64(yaml_obj["if_PN2LN_nACHed"].(int))
	// ...
	scale_slow_PN = yaml_obj["scale_slow_PN"].(float64)
	scale_GABA_PN = yaml_obj["scale_GABA_PN"].(float64)
	scale_GABA_LN = yaml_obj["scale_GABA_LN"].(float64)
	scale_nACH_PN = yaml_obj["scale_nACH_PN"].(float64)
	scale_nACH_LN = yaml_obj["scale_nACH_LN"].(float64)
	// ...
	BG_input_rate = yaml_obj["BG_input_rate"].(float64)
	ORN_input_rate = yaml_obj["ORN_input_rate"].(float64)
	BG_input_strength_PN = yaml_obj["BG_input_strength_PN"].(float64)
	BG_input_strength_LN = yaml_obj["BG_input_strength_LN"].(float64)
	ORN_input_strength_PN = yaml_obj["ORN_input_strength_PN"].(float64)
	ORN_input_strength_LN = yaml_obj["ORN_input_strength_LN"].(float64)
}

func proc_filenames() {
	var randInt int64 = int64(rand.Intn(10e15))
	defer Exit(Enter("$FN()"))
	filename_vRaster_PN = fmt.Sprintf("%s %15d", filename_vRaster_PN, randInt)
	filename_vRaster_LN = fmt.Sprintf("%s %15d", filename_vRaster_LN, randInt)
	filename_stimRaster_PN = fmt.Sprintf("%s %15d", filename_stimRaster_PN, randInt)
	filename_stimRaster_LN = fmt.Sprintf("%s %15d", filename_stimRaster_LN, randInt)
	filename_synaRaster_slow_PN = fmt.Sprintf("%s %15d", filename_synaRaster_slow_PN, randInt)
	filename_synaRaster_GABA_PN = fmt.Sprintf("%s %15d", filename_synaRaster_GABA_PN, randInt)
	filename_synaRaster_GABA_LN = fmt.Sprintf("%s %15d", filename_synaRaster_GABA_LN, randInt)
	filename_synaRaster_nACH_PN = fmt.Sprintf("%s %15d", filename_synaRaster_nACH_PN, randInt)
	filename_synaRaster_nACH_LN = fmt.Sprintf("%s %15d", filename_synaRaster_nACH_LN, randInt)
	// ...
	filename_vFluct_PN = fmt.Sprintf("%s %15d", filename_vFluct_PN, randInt)
	filename_vFluct_LN = fmt.Sprintf("%s %15d", filename_vFluct_LN, randInt)
	filename_stimFluct_PN = fmt.Sprintf("%s %15d", filename_stimFluct_PN, randInt)
	filename_stimFluct_LN = fmt.Sprintf("%s %15d", filename_stimFluct_LN, randInt)
	filename_synaFluct_slow_PN = fmt.Sprintf("%s %15d", filename_synaFluct_slow_PN, randInt)
	filename_synaFluct_GABA_PN = fmt.Sprintf("%s %15d", filename_synaFluct_GABA_PN, randInt)
	filename_synaFluct_GABA_LN = fmt.Sprintf("%s %15d", filename_synaFluct_GABA_LN, randInt)
	filename_synaFluct_nACH_PN = fmt.Sprintf("%s %15d", filename_synaFluct_nACH_PN, randInt)
	filename_synaFluct_nACH_LN = fmt.Sprintf("%s %15d", filename_synaFluct_nACH_LN, randInt)
}

func open_plot_files() {
	defer Exit(Enter("$FN()"))
	file_vRaster_PN, err = os.OpenFile(filename_vRaster_PN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	file_vRaster_LN, err = os.OpenFile(filename_vRaster_LN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	file_stimRaster_PN, err = os.OpenFile(filename_stimRaster_PN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	file_stimRaster_LN, err = os.OpenFile(filename_stimRaster_LN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	// ...
	file_synaRaster_slow_PN, err = os.OpenFile(filename_synaRaster_slow_PN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	file_synaRaster_GABA_PN, err = os.OpenFile(filename_synaRaster_GABA_PN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	file_synaRaster_GABA_LN, err = os.OpenFile(filename_synaRaster_GABA_LN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	file_synaRaster_nACH_PN, err = os.OpenFile(filename_synaRaster_nACH_PN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	file_synaRaster_nACH_LN, err = os.OpenFile(filename_synaRaster_nACH_LN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	// ...
	file_vFluct_PN, err = os.OpenFile(filename_vFluct_PN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	file_vFluct_LN, err = os.OpenFile(filename_vFluct_LN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	file_stimFluct_PN, err = os.OpenFile(filename_stimFluct_PN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	file_stimFluct_LN, err = os.OpenFile(filename_stimFluct_LN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	// ...
	file_synaFluct_slow_PN, err = os.OpenFile(filename_synaFluct_slow_PN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	file_synaFluct_GABA_PN, err = os.OpenFile(filename_synaFluct_GABA_PN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	file_synaFluct_GABA_LN, err = os.OpenFile(filename_synaFluct_GABA_LN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	file_synaFluct_nACH_PN, err = os.OpenFile(filename_synaFluct_nACH_PN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	file_synaFluct_nACH_LN, err = os.OpenFile(filename_synaFluct_nACH_LN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
}

func close_plot_files() {
	defer Exit(Enter("$FN()"))
	file_vFluct_PN.Sync()
	file_vFluct_LN.Sync()
	file_stimFluct_PN.Sync()
	file_stimFluct_LN.Sync()
	// ...
	file_synaFluct_slow_PN.Sync()
	file_synaFluct_GABA_PN.Sync()
	file_synaFluct_GABA_LN.Sync()
	file_synaFluct_nACH_PN.Sync()
	file_synaFluct_nACH_LN.Sync()
	// ...
	file_stimRaster_PN.Sync()
	file_stimRaster_LN.Sync()
	file_vRaster_PN.Sync()
	file_vRaster_LN.Sync()
	// ...
	file_synaRaster_slow_PN.Sync()
	file_synaRaster_GABA_PN.Sync()
	file_synaRaster_GABA_LN.Sync()
	file_synaRaster_nACH_PN.Sync()
	file_synaRaster_nACH_LN.Sync()
	// ...
	file_vFluct_PN.Close()
	file_vFluct_LN.Close()
	file_stimFluct_PN.Close()
	file_stimFluct_LN.Close()
	// ...
	file_synaFluct_slow_PN.Close()
	file_synaFluct_GABA_PN.Close()
	file_synaFluct_GABA_LN.Close()
	file_synaFluct_nACH_PN.Close()
	file_synaFluct_nACH_LN.Close()
	// ...
	file_stimRaster_PN.Close()
	file_stimRaster_LN.Close()
	file_vRaster_PN.Close()
	file_vRaster_LN.Close()
	// ...
	file_synaRaster_slow_PN.Close()
	file_synaRaster_GABA_PN.Close()
	file_synaRaster_GABA_LN.Close()
	file_synaRaster_nACH_PN.Close()
	file_synaRaster_nACH_LN.Close()
}

func init_stimRaster_plot() {
	defer Exit(Enter("$FN()"))
	plotname := ""
	persist := false // if close the plotting frame when the plot is done
	debug := false
	stimPalette := "set palette defined(0 '#FFFFCC', -1 '#FFEDA0', -2 '#FED976', -3 '#FEB24C', -5 '#FD8D3C', -8 '#FC4E2A', -11 '#E31A1C', -14 '#B10026')"
	plot_stimRaster_PN, err = gnuplot.NewPlotter(plotname, persist, debug)
	check(err)
	plot_stimRaster_LN, err = gnuplot.NewPlotter(plotname, persist, debug)
	check(err)
	// ...
	plot_stimRaster_PN.CheckedCmd("reset")
	plot_stimRaster_PN.CheckedCmd("set term gif size 800,400 delay 100")
	plot_stimRaster_PN.CheckedCmd("set term X11")
	plot_stimRaster_PN.CheckedCmd("set title \" Stimulus (uA) feed into each PN\"")
	plot_stimRaster_PN.SetXLabel("Current Time (ms)")
	plot_stimRaster_PN.SetYLabel("PN number")
	plot_stimRaster_PN.CheckedCmd("set cbrange [-2.25:-0.75]")
	plot_stimRaster_PN.CheckedCmd("set view map")
	plot_stimRaster_PN.CheckedCmd(stimPalette)
	plot_stimRaster_PN.CheckedCmd("set yrange [-0.5:%v]", float64(PN_number)-0.5)
	// ...
	plot_stimRaster_LN.CheckedCmd("reset")
	plot_stimRaster_LN.CheckedCmd("set term gif size 800,400 delay 100")
	plot_stimRaster_LN.CheckedCmd("set term X11")
	plot_stimRaster_LN.CheckedCmd("set title \" Stimulus (uA) feed into each LN\"")
	plot_stimRaster_LN.SetXLabel("Current Time (ms)")
	plot_stimRaster_LN.SetYLabel("LN number")
	plot_stimRaster_LN.CheckedCmd("set cbrange [-1:0]")
	plot_stimRaster_LN.CheckedCmd("set view map")
	plot_stimRaster_LN.CheckedCmd(stimPalette)
	plot_stimRaster_LN.CheckedCmd("set yrange [-0.5:%v]", float64(LN_number)-0.5)
}

func init_synaRaster_plot() {
	defer Exit(Enter("$FN()"))
	plotname := ""
	persist := false // if close the plotting frame when the plot is done
	debug := false
	stimPalette := "set palette defined(0 '#FFFFCC', -1 '#FFEDA0', -2 '#FED976', -3 '#FEB24C', -5 '#FD8D3C', -8 '#FC4E2A', -11 '#E31A1C', -14 '#B10026')"
	plot_synaRaster_slow_PN, err = gnuplot.NewPlotter(plotname, persist, debug)
	check(err)
	plot_synaRaster_GABA_PN, err = gnuplot.NewPlotter(plotname, persist, debug)
	check(err)
	plot_synaRaster_GABA_LN, err = gnuplot.NewPlotter(plotname, persist, debug)
	check(err)
	plot_synaRaster_nACH_PN, err = gnuplot.NewPlotter(plotname, persist, debug)
	check(err)
	plot_synaRaster_nACH_LN, err = gnuplot.NewPlotter(plotname, persist, debug)
	check(err)
	// ...
	plot_synaRaster_slow_PN.CheckedCmd("reset")
	plot_synaRaster_slow_PN.CheckedCmd("set term gif size 800,400 delay 100")
	plot_synaRaster_slow_PN.CheckedCmd("set term X11")
	plot_synaRaster_slow_PN.CheckedCmd("set title \" slow synapse current (uA) feed into each PN\"")
	plot_synaRaster_slow_PN.SetXLabel("Current Time (ms)")
	plot_synaRaster_slow_PN.SetYLabel("PN number")
	plot_synaRaster_slow_PN.CheckedCmd("set cbrange [-0.1:2]")
	plot_synaRaster_slow_PN.CheckedCmd("set view map")
	plot_synaRaster_slow_PN.CheckedCmd(stimPalette)
	plot_synaRaster_slow_PN.CheckedCmd("set yrange [-0.5:%v]", float64(PN_number)-0.5)
	// ...
	plot_synaRaster_GABA_PN.CheckedCmd("reset")
	plot_synaRaster_GABA_PN.CheckedCmd("set term gif size 800,400 delay 100")
	plot_synaRaster_GABA_PN.CheckedCmd("set term X11")
	plot_synaRaster_GABA_PN.CheckedCmd("set title \" GABA synapse current (uA) feed into each PN\"")
	plot_synaRaster_GABA_PN.SetXLabel("Current Time (ms)")
	plot_synaRaster_GABA_PN.SetYLabel("PN number")
	plot_synaRaster_GABA_PN.CheckedCmd("set cbrange [-5:20]")
	plot_synaRaster_GABA_PN.CheckedCmd("set view map")
	plot_synaRaster_GABA_PN.CheckedCmd(stimPalette)
	plot_synaRaster_GABA_PN.CheckedCmd("set yrange [-0.5:%v]", float64(PN_number)-0.5)
	// ...
	plot_synaRaster_GABA_LN.CheckedCmd("reset")
	plot_synaRaster_GABA_LN.CheckedCmd("set term gif size 800,400 delay 100")
	plot_synaRaster_GABA_LN.CheckedCmd("set term X11")
	plot_synaRaster_GABA_LN.CheckedCmd("set title \" GABA synapse current (uA) feed into each LN\"")
	plot_synaRaster_GABA_LN.SetXLabel("Current Time (ms)")
	plot_synaRaster_GABA_LN.SetYLabel("LN number")
	plot_synaRaster_GABA_LN.CheckedCmd("set cbrange [-0.1:2]")
	plot_synaRaster_GABA_LN.CheckedCmd("set view map")
	plot_synaRaster_GABA_LN.CheckedCmd(stimPalette)
	plot_synaRaster_GABA_LN.CheckedCmd("set yrange [-0.5:%v]", float64(LN_number)-0.5)
	// ...
	plot_synaRaster_nACH_PN.CheckedCmd("reset")
	plot_synaRaster_nACH_PN.CheckedCmd("set term gif size 800,400 delay 100")
	plot_synaRaster_nACH_PN.CheckedCmd("set term X11")
	plot_synaRaster_nACH_PN.CheckedCmd("set title \" nACH synapse current (uA) feed into each PN\"")
	plot_synaRaster_nACH_PN.SetXLabel("Current Time (ms)")
	plot_synaRaster_nACH_PN.SetYLabel("PN number")
	plot_synaRaster_nACH_PN.CheckedCmd("set cbrange [-2:1]")
	plot_synaRaster_nACH_PN.CheckedCmd("set view map")
	plot_synaRaster_nACH_PN.CheckedCmd(stimPalette)
	plot_synaRaster_nACH_PN.CheckedCmd("set yrange [-0.5:%v]", float64(PN_number)-0.5)
	// ...
	plot_synaRaster_nACH_LN.CheckedCmd("reset")
	plot_synaRaster_nACH_LN.CheckedCmd("set term gif size 800,400 delay 100")
	plot_synaRaster_nACH_LN.CheckedCmd("set term X11")
	plot_synaRaster_nACH_LN.CheckedCmd("set title \" nACH synapse current (uA) feed into each LN\"")
	plot_synaRaster_nACH_LN.SetXLabel("Current Time (ms)")
	plot_synaRaster_nACH_LN.SetYLabel("LN number")
	plot_synaRaster_nACH_LN.CheckedCmd("set cbrange [-10:1]")
	plot_synaRaster_nACH_LN.CheckedCmd("set view map")
	plot_synaRaster_nACH_LN.CheckedCmd(stimPalette)
	plot_synaRaster_nACH_LN.CheckedCmd("set yrange [-0.5:%v]", float64(LN_number)-0.5)
}

func init_vRaster_plot() {
	defer Exit(Enter("$FN()"))
	plotname := ""
	persist := false // if close the plotting frame when the plot is done
	debug := false
	vPalette := "set palette defined(0 '#FFFFCC', 3 '#FFEDA0', 6 '#FED976', 9 '#FEB24C', 11 '#FD8D3C', 12 '#FC4E2A', 13 '#E31A1C', 14 '#B10026')"
	plot_vRaster_PN, err = gnuplot.NewPlotter(plotname, persist, debug)
	check(err)
	plot_vRaster_LN, err = gnuplot.NewPlotter(plotname, persist, debug)
	check(err)
	// ...
	plot_vRaster_PN.CheckedCmd("reset")
	plot_vRaster_PN.CheckedCmd("set term gif size 800,400 delay 100")
	plot_vRaster_PN.CheckedCmd("set term X11")
	plot_vRaster_PN.CheckedCmd("set title \" membrane voltage (mV) raster of each PN\"")
	plot_vRaster_PN.SetXLabel("Current Time (ms)")
	plot_vRaster_PN.SetYLabel("PN number")
	plot_vRaster_PN.CheckedCmd("set cbrange [-80:40]")
	plot_vRaster_PN.CheckedCmd("set view map")
	plot_vRaster_PN.CheckedCmd(vPalette)
	plot_vRaster_PN.CheckedCmd("set yrange [-0.5:%v]", float64(PN_number)-0.5)
	// ...
	plot_vRaster_LN.CheckedCmd("reset")
	plot_vRaster_LN.CheckedCmd("set term gif size 800,400 delay 100")
	plot_vRaster_LN.CheckedCmd("set term X11")
	plot_vRaster_LN.CheckedCmd("set title \" membrane voltage (mV) raster of each LN\"")
	plot_vRaster_LN.SetXLabel("Current Time (ms)")
	plot_vRaster_LN.SetYLabel("LN number")
	plot_vRaster_LN.CheckedCmd("set cbrange [-80:0]")
	plot_vRaster_LN.CheckedCmd("set view map")
	plot_vRaster_LN.CheckedCmd(vPalette)
	plot_vRaster_LN.CheckedCmd("set yrange [-0.5:%v]", float64(LN_number)-0.5)
}

func init_vFluct_plot() {
	defer Exit(Enter("$FN()"))
	plotname := ""
	persist := false // if close the plotting frame when the plot is done
	debug := false
	plot_vFluct_PN, err = gnuplot.NewPlotter(plotname, persist, debug)
	check(err)
	plot_vFluct_LN, err = gnuplot.NewPlotter(plotname, persist, debug)
	check(err)
	// ...
	plot_vFluct_PN.CheckedCmd("reset")
	plot_vFluct_PN.CheckedCmd("set term gif size 800,400 delay 100")
	plot_vFluct_PN.CheckedCmd("set term X11")
	plot_vFluct_PN.CheckedCmd("set title \" membrane voltage (mV) fluction of given PN \"")
	plot_vFluct_PN.SetXLabel("Current Time (ms)")
	plot_vFluct_PN.SetYLabel("Voltage (mV)")
	plot_vFluct_PN.CheckedCmd("set yrange [-80:40]")
	// ...
	plot_vFluct_LN.CheckedCmd("reset")
	plot_vFluct_LN.CheckedCmd("set term gif size 800,400 delay 100")
	plot_vFluct_LN.CheckedCmd("set term X11")
	plot_vFluct_LN.CheckedCmd("set title \" membrane voltage (mV) fluction of given LN \"")
	plot_vFluct_LN.SetXLabel("Current Time (ms)")
	plot_vFluct_LN.SetYLabel("Voltage(mV)")
	plot_vFluct_LN.CheckedCmd("set yrange [-80:0]")
}

func init_stimFluct_plot() {
	defer Exit(Enter("$FN()"))
	plotname := ""
	persist := false // if close the plotting frame when the plot is done
	debug := false
	plot_stimFluct_PN, err = gnuplot.NewPlotter(plotname, persist, debug)
	check(err)
	plot_stimFluct_LN, err = gnuplot.NewPlotter(plotname, persist, debug)
	check(err)
	// ...
	plot_stimFluct_PN.CheckedCmd("reset")
	plot_stimFluct_PN.CheckedCmd("set term gif size 800,400 delay 100")
	plot_stimFluct_PN.CheckedCmd("set term X11")
	plot_stimFluct_PN.CheckedCmd("set title \" stimulus current (uA) fluction to given PN \"")
	plot_stimFluct_PN.SetXLabel("Current Time (ms)")
	plot_stimFluct_PN.SetYLabel("Current (uA)")
	plot_stimFluct_PN.CheckedCmd("set yrange [-3:0]")
	// ...
	plot_stimFluct_LN.CheckedCmd("reset")
	plot_stimFluct_LN.CheckedCmd("set term gif size 800,400 delay 100")
	plot_stimFluct_LN.CheckedCmd("set term X11")
	plot_stimFluct_LN.CheckedCmd("set title \" stimulus current (uA) fluction to given LN \"")
	plot_stimFluct_LN.SetXLabel("Current Time (ms)")
	plot_stimFluct_LN.SetYLabel("Current (uA)")
	plot_stimFluct_LN.CheckedCmd("set yrange [-1:0.1]")
}

func init_synaFluct_plot() {
	defer Exit(Enter("$FN()"))
	plotname := ""
	persist := false // if close the plotting frame when the plot is done
	debug := false
	plot_synaFluct_slow_PN, err = gnuplot.NewPlotter(plotname, persist, debug)
	check(err)
	plot_synaFluct_GABA_PN, err = gnuplot.NewPlotter(plotname, persist, debug)
	check(err)
	plot_synaFluct_GABA_LN, err = gnuplot.NewPlotter(plotname, persist, debug)
	check(err)
	plot_synaFluct_nACH_PN, err = gnuplot.NewPlotter(plotname, persist, debug)
	check(err)
	plot_synaFluct_nACH_LN, err = gnuplot.NewPlotter(plotname, persist, debug)
	check(err)
	// ...
	plot_synaFluct_slow_PN.CheckedCmd("reset")
	plot_synaFluct_slow_PN.CheckedCmd("set term gif size 800,400 delay 100")
	plot_synaFluct_slow_PN.CheckedCmd("set term X11")
	plot_synaFluct_slow_PN.CheckedCmd("set title \" slow Synapse Current (uA) fluction from LNs to given PN \"")
	plot_synaFluct_slow_PN.SetXLabel("Current Time (ms)")
	plot_synaFluct_slow_PN.SetYLabel("Voltage (uA)")
	plot_synaFluct_slow_PN.CheckedCmd("set yrange [-0.1:2]")
	// ...
	plot_synaFluct_GABA_PN.CheckedCmd("reset")
	plot_synaFluct_GABA_PN.CheckedCmd("set term gif size 800,400 delay 100")
	plot_synaFluct_GABA_PN.CheckedCmd("set term X11")
	plot_synaFluct_GABA_PN.CheckedCmd("set title \" GABA Synapse Current (uA) fluction from LNs to given PN \"")
	plot_synaFluct_GABA_PN.SetXLabel("Current Time (ms)")
	plot_synaFluct_GABA_PN.SetYLabel("Voltage (uA)")
	plot_synaFluct_GABA_PN.CheckedCmd("set yrange [-5:20]")
	// ...
	plot_synaFluct_GABA_LN.CheckedCmd("reset")
	plot_synaFluct_GABA_LN.CheckedCmd("set term gif size 800,400 delay 100")
	plot_synaFluct_GABA_LN.CheckedCmd("set term X11")
	plot_synaFluct_GABA_LN.CheckedCmd("set title \" GABA Synapse Current (uA) fluction from LNs to given LN \"")
	plot_synaFluct_GABA_LN.SetXLabel("Current Time (ms)")
	plot_synaFluct_GABA_LN.SetYLabel("Voltage (uA)")
	plot_synaFluct_GABA_LN.CheckedCmd("set yrange [-0.1:2]")
	// ...
	plot_synaFluct_nACH_PN.CheckedCmd("reset")
	plot_synaFluct_nACH_PN.CheckedCmd("set term gif size 800,400 delay 100")
	plot_synaFluct_nACH_PN.CheckedCmd("set term X11")
	plot_synaFluct_nACH_PN.CheckedCmd("set title \" nACH Synapse Current (uA) fluction from PNs to given PN \"")
	plot_synaFluct_nACH_PN.SetXLabel("Current Time (ms)")
	plot_synaFluct_nACH_PN.SetYLabel("Voltage (uA)")
	plot_synaFluct_nACH_PN.CheckedCmd("set yrange [-2:1]")
	// ...
	plot_synaFluct_nACH_LN.CheckedCmd("reset")
	plot_synaFluct_nACH_LN.CheckedCmd("set term gif size 800,400 delay 100")
	plot_synaFluct_nACH_LN.CheckedCmd("set term X11")
	plot_synaFluct_nACH_LN.CheckedCmd("set title \" nACH Synapse Current (uA) fluction from PNs to given LN \"")
	plot_synaFluct_nACH_LN.SetXLabel("Current Time (ms)")
	plot_synaFluct_nACH_LN.SetYLabel("Voltage(uA)")
	plot_synaFluct_nACH_LN.CheckedCmd("set yrange [-10:1]")
}

func init_plot() {
	defer Exit(Enter("$FN()"))
	init_stimRaster_plot()
	init_vRaster_plot()
	init_synaRaster_plot()
	init_vFluct_plot()
	init_stimFluct_plot()
	init_synaFluct_plot()
}

func dest_plot() {
	defer Exit(Enter("$FN()"))
	plot_stimRaster_PN.CheckedCmd("set output")
	plot_stimRaster_PN.CheckedCmd("replot")
	plot_stimRaster_PN.CheckedCmd("quit")
	if exist_file(filename_stimRaster_PN) {
		os.Remove(filename_stimRaster_PN)
	}
	plot_stimRaster_LN.CheckedCmd("set output")
	plot_stimRaster_LN.CheckedCmd("replot")
	plot_stimRaster_LN.CheckedCmd("quit")
	if exist_file(filename_stimRaster_LN) {
		os.Remove(filename_stimRaster_LN)
	}
	// ...
	plot_vRaster_PN.CheckedCmd("set output")
	plot_vRaster_PN.CheckedCmd("replot")
	plot_vRaster_PN.CheckedCmd("quit")
	if exist_file(filename_vRaster_PN) {
		os.Remove(filename_vRaster_PN)
	}
	plot_vRaster_LN.CheckedCmd("set output")
	plot_vRaster_LN.CheckedCmd("replot")
	plot_vRaster_LN.CheckedCmd("quit")
	if exist_file(filename_vRaster_LN) {
		os.Remove(filename_vRaster_LN)
	}
	// ...
	plot_synaRaster_slow_PN.CheckedCmd("set output")
	plot_synaRaster_slow_PN.CheckedCmd("replot")
	plot_synaRaster_slow_PN.CheckedCmd("quit")
	if exist_file(filename_synaRaster_slow_PN) {
		os.Remove(filename_synaRaster_slow_PN)
	}
	plot_synaRaster_GABA_PN.CheckedCmd("set output")
	plot_synaRaster_GABA_PN.CheckedCmd("replot")
	plot_synaRaster_GABA_PN.CheckedCmd("quit")
	if exist_file(filename_synaRaster_GABA_PN) {
		os.Remove(filename_synaRaster_GABA_PN)
	}
	plot_synaRaster_GABA_LN.CheckedCmd("set output")
	plot_synaRaster_GABA_LN.CheckedCmd("replot")
	plot_synaRaster_GABA_LN.CheckedCmd("quit")
	if exist_file(filename_synaRaster_GABA_LN) {
		os.Remove(filename_synaRaster_GABA_LN)
	}
	plot_synaRaster_nACH_PN.CheckedCmd("set output")
	plot_synaRaster_nACH_PN.CheckedCmd("replot")
	plot_synaRaster_nACH_PN.CheckedCmd("quit")
	if exist_file(filename_synaRaster_nACH_PN) {
		os.Remove(filename_synaRaster_nACH_PN)
	}
	plot_synaRaster_nACH_LN.CheckedCmd("set output")
	plot_synaRaster_nACH_LN.CheckedCmd("replot")
	plot_synaRaster_nACH_LN.CheckedCmd("quit")
	if exist_file(filename_synaRaster_nACH_LN) {
		os.Remove(filename_synaRaster_nACH_LN)
	}
	// ...
	plot_vFluct_PN.CheckedCmd("set output")
	plot_vFluct_PN.CheckedCmd("replot")
	plot_vFluct_PN.CheckedCmd("quit")
	if exist_file(filename_vFluct_PN) {
		os.Remove(filename_vFluct_PN)
	}
	plot_vFluct_LN.CheckedCmd("set output")
	plot_vFluct_LN.CheckedCmd("replot")
	plot_vFluct_LN.CheckedCmd("quit")
	if exist_file(filename_vFluct_LN) {
		os.Remove(filename_vFluct_LN)
	}
	// ...
	plot_stimFluct_PN.CheckedCmd("set output")
	plot_stimFluct_PN.CheckedCmd("replot")
	plot_stimFluct_PN.CheckedCmd("quit")
	if exist_file(filename_stimFluct_PN) {
		os.Remove(filename_stimFluct_PN)
	}
	plot_stimFluct_LN.CheckedCmd("set output")
	plot_stimFluct_LN.CheckedCmd("replot")
	plot_stimFluct_LN.CheckedCmd("quit")
	if exist_file(filename_stimFluct_LN) {
		os.Remove(filename_stimFluct_LN)
	}
	// ...
	plot_synaFluct_slow_PN.CheckedCmd("set output")
	plot_synaFluct_slow_PN.CheckedCmd("replot")
	plot_synaFluct_slow_PN.CheckedCmd("quit")
	if exist_file(filename_synaFluct_slow_PN) {
		os.Remove(filename_synaFluct_slow_PN)
	}
	plot_synaFluct_GABA_PN.CheckedCmd("set output")
	plot_synaFluct_GABA_PN.CheckedCmd("replot")
	plot_synaFluct_GABA_PN.CheckedCmd("quit")
	if exist_file(filename_synaFluct_GABA_PN) {
		os.Remove(filename_synaFluct_GABA_PN)
	}
	plot_synaFluct_GABA_LN.CheckedCmd("set output")
	plot_synaFluct_GABA_LN.CheckedCmd("replot")
	plot_synaFluct_GABA_LN.CheckedCmd("quit")
	if exist_file(filename_synaFluct_GABA_LN) {
		os.Remove(filename_synaFluct_GABA_LN)
	}
	plot_synaFluct_nACH_PN.CheckedCmd("set output")
	plot_synaFluct_nACH_PN.CheckedCmd("replot")
	plot_synaFluct_nACH_PN.CheckedCmd("quit")
	if exist_file(filename_synaFluct_nACH_PN) {
		os.Remove(filename_synaFluct_nACH_PN)
	}
	plot_synaFluct_nACH_LN.CheckedCmd("set output")
	plot_synaFluct_nACH_LN.CheckedCmd("replot")
	plot_synaFluct_nACH_LN.CheckedCmd("quit")
	if exist_file(filename_synaFluct_nACH_LN) {
		os.Remove(filename_synaFluct_nACH_LN)
	}
}

func write_plot_files_PN(cur_ms int64) { // update the data that are need to written and plot
	var id, plotStep int64
	defer Exit(Enter("$FN(%v)", cur_ms))
	if if_PN_plot <= 0 {
		return
	}
	for plotStep = int64(math.Max(0.0, float64(cur_ms-ms_per_frame))); plotStep < cur_ms; plotStep += 1 {
		if if_stimRaster_plot > 0 {
			for id = 0; id < PN_number; id++ {
				_, _ = file_stimRaster_PN.Write([]byte(strconv.FormatFloat(float64(plotStep), 'f', 15, 64)))
				_, _ = file_stimRaster_PN.WriteString("  ")
				_, _ = file_stimRaster_PN.Write([]byte(strconv.FormatInt(id, 10))) // 10 based
				_, _ = file_stimRaster_PN.WriteString("  ")
				_, _ = file_stimRaster_PN.Write([]byte(strconv.FormatFloat(PN_stim_list[plotStep][id], 'f', 15, 64))) // 15 digits, 64 bits
				_, _ = file_stimRaster_PN.WriteString("\n")
			}
			_, _ = file_stimRaster_PN.WriteString("\n")
		}
		// write both stim and v files...
		if if_vRaster_plot > 0 {
			for id = 0; id < PN_number; id++ {
				_, _ = file_vRaster_PN.Write([]byte(strconv.FormatFloat(float64(plotStep), 'f', 15, 64)))
				_, _ = file_vRaster_PN.WriteString("  ")
				_, _ = file_vRaster_PN.Write([]byte(strconv.FormatInt(id, 10)))
				_, _ = file_vRaster_PN.WriteString("  ")
				_, _ = file_vRaster_PN.Write([]byte(strconv.FormatFloat(PN_V_list[plotStep][id], 'f', 15, 64)))
				_, _ = file_vRaster_PN.WriteString("\n")
			}
			_, _ = file_vRaster_PN.WriteString("\n")
		}
		// ...
		if if_synaRaster_plot > 0 {
			for id = 0; id < PN_number; id++ {
				_, _ = file_synaRaster_slow_PN.Write([]byte(strconv.FormatFloat(float64(plotStep), 'f', 15, 64)))
				_, _ = file_synaRaster_slow_PN.WriteString("  ")
				_, _ = file_synaRaster_slow_PN.Write([]byte(strconv.FormatInt(id, 10))) // 10 based
				_, _ = file_synaRaster_slow_PN.WriteString("  ")
				_, _ = file_synaRaster_slow_PN.Write([]byte(strconv.FormatFloat(PN_stim_list[plotStep][id], 'f', 15, 64)))
				_, _ = file_synaRaster_slow_PN.WriteString("\n")
				// ...
				_, _ = file_synaRaster_GABA_PN.Write([]byte(strconv.FormatFloat(float64(plotStep), 'f', 15, 64)))
				_, _ = file_synaRaster_GABA_PN.WriteString("  ")
				_, _ = file_synaRaster_GABA_PN.Write([]byte(strconv.FormatInt(id, 10))) // 10 based
				_, _ = file_synaRaster_GABA_PN.WriteString("  ")
				_, _ = file_synaRaster_GABA_PN.Write([]byte(strconv.FormatFloat(PN_stim_list[plotStep][id], 'f', 15, 64)))
				_, _ = file_synaRaster_GABA_PN.WriteString("\n")
				// ...
				_, _ = file_synaRaster_nACH_PN.Write([]byte(strconv.FormatFloat(float64(plotStep), 'f', 15, 64)))
				_, _ = file_synaRaster_nACH_PN.WriteString("  ")
				_, _ = file_synaRaster_nACH_PN.Write([]byte(strconv.FormatInt(id, 10))) // 10 based
				_, _ = file_synaRaster_nACH_PN.WriteString("  ")
				_, _ = file_synaRaster_nACH_PN.Write([]byte(strconv.FormatFloat(PN_stim_list[plotStep][id], 'f', 15, 64)))
				_, _ = file_synaRaster_nACH_PN.WriteString("\n")
			}
			_, _ = file_synaRaster_slow_PN.WriteString("\n")
			_, _ = file_synaRaster_GABA_PN.WriteString("\n")
			_, _ = file_synaRaster_nACH_PN.WriteString("\n")
		}
	}
}

func write_plot_files_LN(cur_ms int64) { // update the data that are need to written and plot
	var id, plotStep int64
	defer Exit(Enter("$FN(%v)", cur_ms))
	if if_LN_plot <= 0 {
		return
	}
	for plotStep = int64(math.Max(0.0, float64(cur_ms-ms_per_frame))); plotStep < cur_ms; plotStep += 1 {
		if if_stimRaster_plot > 0 {
			for id = 0; id < LN_number; id++ {
				_, _ = file_stimRaster_LN.Write([]byte(strconv.FormatFloat(float64(plotStep), 'f', 15, 64)))
				_, _ = file_stimRaster_LN.WriteString("  ")
				_, _ = file_stimRaster_LN.Write([]byte(strconv.FormatInt(id, 10))) // 10 based
				_, _ = file_stimRaster_LN.WriteString("  ")
				_, _ = file_stimRaster_LN.Write([]byte(strconv.FormatFloat(LN_stim_list[plotStep][id], 'f', 15, 64))) // 15 digits, 64 bits
				_, _ = file_stimRaster_LN.WriteString("\n")
			}
			_, _ = file_stimRaster_LN.WriteString("\n")
		}
		// write both stim and v files...
		if if_vRaster_plot > 0 {
			for id = 0; id < LN_number; id++ {
				_, _ = file_vRaster_LN.Write([]byte(strconv.FormatFloat(float64(plotStep), 'f', 15, 64)))
				_, _ = file_vRaster_LN.WriteString("  ")
				_, _ = file_vRaster_LN.Write([]byte(strconv.FormatInt(id, 10)))
				_, _ = file_vRaster_LN.WriteString("  ")
				_, _ = file_vRaster_LN.Write([]byte(strconv.FormatFloat(LN_V_list[plotStep][id], 'f', 15, 64)))
				_, _ = file_vRaster_LN.WriteString("\n")
			}
			_, _ = file_vRaster_LN.WriteString("\n")
		}
		// ...
		if if_synaRaster_plot > 0 {
			for id = 0; id < LN_number; id++ {
				_, _ = file_synaRaster_GABA_LN.Write([]byte(strconv.FormatFloat(float64(plotStep), 'f', 15, 64)))
				_, _ = file_synaRaster_GABA_LN.WriteString("  ")
				_, _ = file_synaRaster_GABA_LN.Write([]byte(strconv.FormatInt(id, 10))) // 10 based
				_, _ = file_synaRaster_GABA_LN.WriteString("  ")
				_, _ = file_synaRaster_GABA_LN.Write([]byte(strconv.FormatFloat(LN_stim_list[plotStep][id], 'f', 15, 64)))
				_, _ = file_synaRaster_GABA_LN.WriteString("\n")
				// ...
				_, _ = file_synaRaster_nACH_LN.Write([]byte(strconv.FormatFloat(float64(plotStep), 'f', 15, 64)))
				_, _ = file_synaRaster_nACH_LN.WriteString("  ")
				_, _ = file_synaRaster_nACH_LN.Write([]byte(strconv.FormatInt(id, 10))) // 10 based
				_, _ = file_synaRaster_nACH_LN.WriteString("  ")
				_, _ = file_synaRaster_nACH_LN.Write([]byte(strconv.FormatFloat(LN_stim_list[plotStep][id], 'f', 15, 64)))
				_, _ = file_synaRaster_nACH_LN.WriteString("\n")
			}
			_, _ = file_synaRaster_GABA_LN.WriteString("\n")
			_, _ = file_synaRaster_nACH_LN.WriteString("\n")
		}
	}
}

func write_vFluct_files(cur_ms int64) { // synchronized with V_xN and stim_xN
	var id, plotStep int64
	defer Exit(Enter("$FN(%v)", cur_ms))
	for plotStep = int64(math.Max(0.0, float64(cur_ms-ms_per_frame))); plotStep < cur_ms; plotStep += 1 {
		if if_vFluct_plot > 0 && if_PN_plot > 0 {
			id = plot_fluct_PN_id
			_, _ = file_vFluct_PN.Write([]byte(strconv.FormatFloat(float64(plotStep), 'f', 15, 64)))
			_, _ = file_vFluct_PN.WriteString("  ")
			_, _ = file_vFluct_PN.Write([]byte(strconv.FormatFloat(PN_V_list[plotStep][id], 'f', 15, 64)))
			_, _ = file_vFluct_PN.WriteString("\n")
		}
		// write both stim and v files...
		if if_vFluct_plot > 0 && if_LN_plot > 0 {
			id = plot_fluct_LN_id
			_, _ = file_vFluct_LN.Write([]byte(strconv.FormatFloat(float64(plotStep), 'f', 15, 64)))
			_, _ = file_vFluct_LN.WriteString("  ")
			_, _ = file_vFluct_LN.Write([]byte(strconv.FormatFloat(LN_V_list[plotStep][id], 'f', 15, 64)))
			_, _ = file_vFluct_LN.WriteString("\n")
		}
	}
}

func write_stimFluct_files(cur_ms int64) { // synchronized with V_xN and stim_xN
	var id, plotStep int64
	defer Exit(Enter("$FN(%v)", cur_ms))
	for plotStep = int64(math.Max(0.0, float64(cur_ms-ms_per_frame))); plotStep < cur_ms; plotStep += 1 {
		if if_stimFluct_plot > 0 && if_PN_plot > 0 {
			id = plot_fluct_PN_id
			_, _ = file_stimFluct_PN.Write([]byte(strconv.FormatFloat(float64(plotStep), 'f', 15, 64)))
			_, _ = file_stimFluct_PN.WriteString("  ")
			_, _ = file_stimFluct_PN.Write([]byte(strconv.FormatFloat(PN_stim_list[plotStep][id], 'f', 15, 64)))
			_, _ = file_stimFluct_PN.WriteString("\n")
		}
		// write both stim and v files...
		if if_stimFluct_plot > 0 && if_LN_plot > 0 {
			id = plot_fluct_LN_id
			_, _ = file_stimFluct_LN.Write([]byte(strconv.FormatFloat(float64(plotStep), 'f', 15, 64)))
			_, _ = file_stimFluct_LN.WriteString("  ")
			_, _ = file_stimFluct_LN.Write([]byte(strconv.FormatFloat(LN_stim_list[plotStep][id], 'f', 15, 64)))
			_, _ = file_stimFluct_LN.WriteString("\n")
		}
	}
}

func write_synaFluct_files(cur_ms int64) { // synchronized with V_xN and stim_xN
	var id, plotStep int64
	defer Exit(Enter("$FN(%v)", cur_ms))
	for plotStep = int64(math.Max(0.0, float64(cur_ms-ms_per_frame))); plotStep < cur_ms; plotStep += 1 {
		if if_synaFluct_plot > 0 && if_PN_plot > 0 { // plot synapses goto PNs
			id = plot_fluct_PN_id
			_, _ = file_synaFluct_slow_PN.Write([]byte(strconv.FormatFloat(float64(plotStep), 'f', 15, 64)))
			_, _ = file_synaFluct_slow_PN.WriteString("  ")
			_, _ = file_synaFluct_slow_PN.Write([]byte(strconv.FormatFloat(PN_slow_list[plotStep][id], 'f', 15, 64)))
			_, _ = file_synaFluct_slow_PN.WriteString("\n")
			// ...
			id = plot_fluct_PN_id
			_, _ = file_synaFluct_GABA_PN.Write([]byte(strconv.FormatFloat(float64(plotStep), 'f', 15, 64)))
			_, _ = file_synaFluct_GABA_PN.WriteString("  ")
			_, _ = file_synaFluct_GABA_PN.Write([]byte(strconv.FormatFloat(PN_GABA_list[plotStep][id], 'f', 15, 64)))
			_, _ = file_synaFluct_GABA_PN.WriteString("\n")
			// ...
			id = plot_fluct_PN_id
			_, _ = file_synaFluct_nACH_PN.Write([]byte(strconv.FormatFloat(float64(plotStep), 'f', 15, 64)))
			_, _ = file_synaFluct_nACH_PN.WriteString("  ")
			_, _ = file_synaFluct_nACH_PN.Write([]byte(strconv.FormatFloat(PN_nACH_list[plotStep][id], 'f', 15, 64)))
			_, _ = file_synaFluct_nACH_PN.WriteString("\n")
		}
		// ...
		if if_synaFluct_plot > 0 && if_LN_plot > 0 { // plot synapses goto LNs
			id = plot_fluct_LN_id
			_, _ = file_synaFluct_GABA_LN.Write([]byte(strconv.FormatFloat(float64(plotStep), 'f', 15, 64)))
			_, _ = file_synaFluct_GABA_LN.WriteString("  ")
			_, _ = file_synaFluct_GABA_LN.Write([]byte(strconv.FormatFloat(LN_GABA_list[plotStep][id], 'f', 15, 64)))
			_, _ = file_synaFluct_GABA_LN.WriteString("\n")
			// ...
			id = plot_fluct_LN_id
			_, _ = file_synaFluct_nACH_LN.Write([]byte(strconv.FormatFloat(float64(plotStep), 'f', 15, 64)))
			_, _ = file_synaFluct_nACH_LN.WriteString("  ")
			_, _ = file_synaFluct_nACH_LN.Write([]byte(strconv.FormatFloat(LN_nACH_list[plotStep][id], 'f', 15, 64)))
			_, _ = file_synaFluct_nACH_LN.WriteString("\n")
		}
	}
}

func proc_plot_files(click int64) {
	defer Exit(Enter("$FN(%v)", click))
	open_plot_files()
	write_plot_files_PN(click)   // write all the points that need to be redrawed
	write_plot_files_LN(click)   // write all the points that need to be redrawed
	write_vFluct_files(click)    // write all the points that need to be redrawed
	write_stimFluct_files(click) // write all the points that need to be redrawed
	write_synaFluct_files(click) // write all the points that need to be redrawed
	close_plot_files()
}

func do_plot(cur_ms int64) {
	defer Exit(Enter("$FN(%v)", cur_ms))
	if if_stimRaster_plot > 0 {
		if if_PN_plot > 0 {
			plot_stimRaster_PN.CheckedCmd("set xrange [%v:%v]", cur_ms-ms_per_frame, cur_ms)
			plot_stimRaster_PN.CheckedCmd("splot \"%v\" u 1:2:3 w p ps 0.5 pt 7 palette notitle", filename_stimRaster_PN)
		}
		if if_LN_plot > 0 {
			plot_stimRaster_LN.CheckedCmd("set xrange [%v:%v]", cur_ms-ms_per_frame, cur_ms)
			plot_stimRaster_LN.CheckedCmd("splot \"%v\" u 1:2:3 w p ps 0.5 pt 7 palette notitle", filename_stimRaster_LN)
		}
	}
	// ...
	if if_vRaster_plot > 0 {
		if if_PN_plot > 0 {
			plot_vRaster_PN.CheckedCmd("set xrange [%v:%v]", cur_ms-ms_per_frame, cur_ms)
			plot_vRaster_PN.CheckedCmd("splot \"%v\" u 1:2:3 w p ps 0.5 pt 7 palette notitle", filename_vRaster_PN)
		}
		if if_LN_plot > 0 {
			plot_vRaster_LN.CheckedCmd("set xrange [%v:%v]", cur_ms-ms_per_frame, cur_ms)
			plot_vRaster_LN.CheckedCmd("splot \"%v\" u 1:2:3 w p ps 0.5 pt 7 palette notitle", filename_vRaster_LN)
		}
	}
	// ...
	if if_synaRaster_plot > 0 {
		if if_PN_plot > 0 {
			plot_synaRaster_slow_PN.CheckedCmd("set xrange [%v:%v]", cur_ms-ms_per_frame, cur_ms)
			plot_synaRaster_slow_PN.CheckedCmd("splot \"%v\" u 1:2:3 w p ps 0.5 pt 7 palette notitle", filename_synaRaster_slow_PN)
			plot_synaRaster_GABA_PN.CheckedCmd("set xrange [%v:%v]", cur_ms-ms_per_frame, cur_ms)
			plot_synaRaster_GABA_PN.CheckedCmd("splot \"%v\" u 1:2:3 w p ps 0.5 pt 7 palette notitle", filename_synaRaster_GABA_PN)
			plot_synaRaster_nACH_PN.CheckedCmd("set xrange [%v:%v]", cur_ms-ms_per_frame, cur_ms)
			plot_synaRaster_nACH_PN.CheckedCmd("splot \"%v\" u 1:2:3 w p ps 0.5 pt 7 palette notitle", filename_synaRaster_nACH_PN)
		}
		if if_LN_plot > 0 {
			plot_synaRaster_GABA_LN.CheckedCmd("set xrange [%v:%v]", cur_ms-ms_per_frame, cur_ms)
			plot_synaRaster_GABA_LN.CheckedCmd("splot \"%v\" u 1:2:3 w p ps 0.5 pt 7 palette notitle", filename_synaRaster_GABA_LN)
			plot_synaRaster_nACH_LN.CheckedCmd("set xrange [%v:%v]", cur_ms-ms_per_frame, cur_ms)
			plot_synaRaster_nACH_LN.CheckedCmd("splot \"%v\" u 1:2:3 w p ps 0.5 pt 7 palette notitle", filename_synaRaster_nACH_LN)
		}
	}
	// ...
	if if_vFluct_plot > 0 {
		if if_PN_plot > 0 {
			plot_vFluct_PN.CheckedCmd("set xrange [%v:%v]", cur_ms-ms_per_frame, cur_ms)
			plot_vFluct_PN.CheckedCmd("plot \"%v\" w l notitle", filename_vFluct_PN)
		}
		if if_LN_plot > 0 {
			plot_vFluct_LN.CheckedCmd("set xrange [%v:%v]", cur_ms-ms_per_frame, cur_ms)
			plot_vFluct_LN.CheckedCmd("plot \"%v\" w l notitle", filename_vFluct_LN)
		}
	}
	// ...
	if if_stimFluct_plot > 0 {
		if if_PN_plot > 0 {
			plot_stimFluct_PN.CheckedCmd("set xrange [%v:%v]", cur_ms-ms_per_frame, cur_ms)
			plot_stimFluct_PN.CheckedCmd("plot \"%v\" w l notitle", filename_stimFluct_PN)
		}
		if if_LN_plot > 0 {
			plot_stimFluct_LN.CheckedCmd("set xrange [%v:%v]", cur_ms-ms_per_frame, cur_ms)
			plot_stimFluct_LN.CheckedCmd("plot \"%v\" w l notitle", filename_stimFluct_LN)
		}
	}
	// ...
	if if_synaFluct_plot > 0 {
		if if_PN_plot > 0 { // synapses goto PNs
			plot_synaFluct_slow_PN.CheckedCmd("set xrange [%v:%v]", cur_ms-ms_per_frame, cur_ms)
			plot_synaFluct_slow_PN.CheckedCmd("plot \"%v\" w l notitle", filename_synaFluct_slow_PN)
			plot_synaFluct_GABA_PN.CheckedCmd("set xrange [%v:%v]", cur_ms-ms_per_frame, cur_ms)
			plot_synaFluct_GABA_PN.CheckedCmd("plot \"%v\" w l notitle", filename_synaFluct_GABA_PN)
			plot_synaFluct_nACH_PN.CheckedCmd("set xrange [%v:%v]", cur_ms-ms_per_frame, cur_ms)
			plot_synaFluct_nACH_PN.CheckedCmd("plot \"%v\" w l notitle", filename_synaFluct_nACH_PN)
		}
		if if_LN_plot > 0 { // synapses goto LNs
			plot_synaFluct_GABA_LN.CheckedCmd("set xrange [%v:%v]", cur_ms-ms_per_frame, cur_ms)
			plot_synaFluct_GABA_LN.CheckedCmd("plot \"%v\" w l notitle", filename_synaFluct_GABA_LN)
			plot_synaFluct_nACH_LN.CheckedCmd("set xrange [%v:%v]", cur_ms-ms_per_frame, cur_ms)
			plot_synaFluct_nACH_LN.CheckedCmd("plot \"%v\" w l notitle", filename_synaFluct_nACH_LN)
		}
	}
}

func init_neuron() {
	var id int64
	defer Exit(Enter("$FN()"))
	for id = 0; id < PN_number; id++ {
		PN_reactor[id][0].V = V_rest_PN + 10*rand.Float64() - 5.0
		PN_reactor[id][0].stim = -1.5
		PN_reactor[id][0].h_K = 1
		PN_reactor[id][0].h_Na = 1
		PN_reactor[id][0].h_A = 1
		PN_reactor[id][0].lastSpike = -10
		PN_V_list[0][id] = PN_reactor[id][0].V
		PN_stim_list[0][id] = PN_reactor[id][0].stim
	} // end for PN
	for id = 0; id < LN_number; id++ {
		LN_reactor[id][0].V = V_rest_LN + 10*rand.Float64() - 5.0
		LN_reactor[id][0].stim = -0.05
		LN_reactor[id][0].h_K = 1
		LN_reactor[id][0].h_Ca = 1
		LN_reactor[id][0].h_B = 1
		LN_reactor[id][0].R_slow = 0.3 * rand.Float64()
		LN_reactor[id][0].lastSpike = -10
		LN_V_list[0][id] = LN_reactor[id][0].V
		LN_stim_list[0][id] = LN_reactor[id][0].stim
	} // end for LN
	// currents are not necessary to be initialized...
}

func save_docs() {
	var i, id int64
	defer Exit(Enter("$FN()"))
	doc_V_PN, err = os.OpenFile(docname_V_PN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	doc_V_LN, err = os.OpenFile(docname_V_LN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	doc_stim_PN, err = os.OpenFile(docname_stim_PN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	doc_stim_LN, err = os.OpenFile(docname_stim_LN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	// ...
	doc_slow_PN, err = os.OpenFile(docname_slow_PN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	doc_GABA_PN, err = os.OpenFile(docname_GABA_PN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	doc_GABA_LN, err = os.OpenFile(docname_GABA_LN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	doc_nACH_PN, err = os.OpenFile(docname_nACH_PN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	doc_nACH_LN, err = os.OpenFile(docname_nACH_LN, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	// ...
	for i = 0; i <= time_end; i++ {
		_, _ = doc_V_PN.Write([]byte(strconv.FormatInt(i, 10))) // current time (ms)
		for id = 0; id < PN_number; id++ {
			_, _ = doc_V_PN.WriteString(" ")
			_, _ = doc_V_PN.Write([]byte(strconv.FormatFloat(PN_V_list[i][id], 'f', 15, 64)))
		}
		_, _ = doc_V_PN.WriteString("\n")
		// ...
		_, _ = doc_stim_PN.Write([]byte(strconv.FormatInt(i, 10)))
		for id = 0; id < PN_number; id++ {
			_, _ = doc_stim_PN.WriteString(" ")
			_, _ = doc_stim_PN.Write([]byte(strconv.FormatFloat(PN_stim_list[i][id], 'f', 15, 64)))
		}
		_, _ = doc_stim_PN.WriteString("\n")
		// ...
		_, _ = doc_V_LN.Write([]byte(strconv.FormatInt(i, 10)))
		for id = 0; id < LN_number; id++ {
			_, _ = doc_V_LN.WriteString(" ")
			_, _ = doc_V_LN.Write([]byte(strconv.FormatFloat(LN_V_list[i][id], 'f', 15, 64)))
		}
		_, _ = doc_V_LN.WriteString("\n")
		// ...
		_, _ = doc_stim_LN.Write([]byte(strconv.FormatInt(i, 10)))
		for id = 0; id < LN_number; id++ {
			_, _ = doc_stim_LN.WriteString(" ")
			_, _ = doc_stim_LN.Write([]byte(strconv.FormatFloat(LN_stim_list[i][id], 'f', 15, 64)))
		}
		_, _ = doc_stim_LN.WriteString("\n")
		// ...
		_, _ = doc_slow_PN.Write([]byte(strconv.FormatInt(i, 10))) // current time (ms)
		for id = 0; id < PN_number; id++ {
			_, _ = doc_slow_PN.WriteString(" ")
			_, _ = doc_slow_PN.Write([]byte(strconv.FormatFloat(PN_slow_list[i][id], 'f', 15, 64)))
		}
		_, _ = doc_slow_PN.WriteString("\n")
		// ...
		_, _ = doc_GABA_PN.Write([]byte(strconv.FormatInt(i, 10))) // current time (ms)
		for id = 0; id < PN_number; id++ {
			_, _ = doc_GABA_PN.WriteString(" ")
			_, _ = doc_GABA_PN.Write([]byte(strconv.FormatFloat(PN_GABA_list[i][id], 'f', 15, 64)))
		}
		_, _ = doc_GABA_PN.WriteString("\n")
		// ...
		_, _ = doc_GABA_LN.Write([]byte(strconv.FormatInt(i, 10)))
		for id = 0; id < LN_number; id++ {
			_, _ = doc_GABA_LN.WriteString(" ")
			_, _ = doc_GABA_LN.Write([]byte(strconv.FormatFloat(LN_GABA_list[i][id], 'f', 15, 64)))
		}
		_, _ = doc_GABA_LN.WriteString("\n")
		// ...
		_, _ = doc_nACH_PN.Write([]byte(strconv.FormatInt(i, 10)))
		for id = 0; id < PN_number; id++ {
			_, _ = doc_nACH_PN.WriteString(" ")
			_, _ = doc_nACH_PN.Write([]byte(strconv.FormatFloat(PN_nACH_list[i][id], 'f', 15, 64)))
		}
		_, _ = doc_nACH_PN.WriteString("\n")
		// ...
		_, _ = doc_nACH_LN.Write([]byte(strconv.FormatInt(i, 10)))
		for id = 0; id < LN_number; id++ {
			_, _ = doc_nACH_LN.WriteString(" ")
			_, _ = doc_nACH_LN.Write([]byte(strconv.FormatFloat(LN_nACH_list[i][id], 'f', 15, 64)))
		}
		_, _ = doc_nACH_LN.WriteString("\n")
	}
	// ...
	doc_stim_PN.Close()
	doc_stim_LN.Close()
	doc_V_PN.Close()
	doc_V_LN.Close()
	doc_slow_PN.Close()
	doc_GABA_PN.Close()
	doc_GABA_LN.Close()
	doc_nACH_PN.Close()
	doc_nACH_LN.Close()
}

func usage() {
	defer Exit(Enter("$FN()"))
	fmt.Println("Usage: awesome [odor_slip] [matrix_seed]")
	fmt.Println("    where odor_slip and matrix_seed are integers.")
	fmt.Println()
}

func scale_synapse_currents() {
	defer Exit(Enter("$FN()"))
	g_GABA_PN *= scale_GABA_PN // LN --[GABA]-> PN --- use this one
	g_nACH_PN *= scale_nACH_PN // PN --[nACH]-> PN
	g_slow_PN *= scale_slow_PN // LN --[slow]-> PN --- use this one
	g_GABA_LN *= scale_GABA_LN // LN --[GABA]-> LN
	g_nACH_LN *= scale_nACH_LN // PN --[nACH]-> LN
	// ...
	fmt.Println("Scaling:")
	fmt.Println("scale_GABA_PN ", scale_GABA_PN)
	fmt.Println("scale_nACH_PN ", scale_nACH_PN)
	fmt.Println("scale_slow_PN ", scale_slow_PN)
	fmt.Println("scale_GABA_LN ", scale_GABA_LN)
	fmt.Println("scale_nACH_LN ", scale_nACH_LN)
	fmt.Println()
}

func remove_coupling_matrix() {
	defer Exit(Enter("$FN()"))
	var i, j int64
	if if_PN2PN_nACHed <= 0 {
		for i = 0; i < PN_number; i++ { // pn to pn
			for j = 0; j < PN_number; j++ {
				PN2PN_nACH_mat[i][j] = 0
			}
		}
	}
	// ...
	if if_LN2PN_GABAed <= 0 {
		for i = 0; i < LN_number; i++ { // ln to pn
			for j = 0; j < PN_number; j++ {
				LN2PN_GABA_mat[i][j] = 0
			}
		}
	}
	// ...
	if if_LN2PN_slowed <= 0 {
		for i = 0; i < LN_number; i++ { // ln to pn
			for j = 0; j < PN_number; j++ {
				LN2PN_slow_mat[i][j] = 0
			}
		}
	}
	// ...
	if if_LN2LN_GABAed <= 0 {
		for i = 0; i < LN_number; i++ { // ln to ln
			for j = 0; j < LN_number; j++ {
				LN2LN_GABA_mat[i][j] = 0
			}
		}
	}
	// ...
	if if_PN2LN_nACHed <= 0 {
		for i = 0; i < PN_number; i++ { // pn to ln
			for j = 0; j < LN_number; j++ {
				PN2LN_nACH_mat[i][j] = 0
			}
		}
	}
	// ...
	fmt.Println("Decouplings:")
	fmt.Println("if_PN2PN_nACHed ", if_PN2PN_nACHed > 0)
	fmt.Println("if_PN2LN_nACHed ", if_PN2LN_nACHed > 0)
	fmt.Println("if_LN2PN_slowed ", if_LN2PN_slowed > 0)
	fmt.Println("if_LN2PN_GABAed ", if_LN2PN_GABAed > 0)
	fmt.Println("if_LN2LN_GABAed ", if_LN2LN_GABAed > 0)
	fmt.Println()
}

func disp_para() {
	defer Exit(Enter("$FN()"))
	fmt.Println("matrix_seed ", matrix_seed)
	fmt.Println("PN_number ", PN_number)
	fmt.Println("LN_number ", LN_number)
	fmt.Println("stim_PN_num ", stim_PN_num)
	fmt.Println("stim_LN_num ", stim_LN_num)
	fmt.Println("ORN_number ", ORN_number)
	fmt.Println("stim_onset ", stim_onset/float64(ms_per_second))
	fmt.Println("stim_rise ", stim_rise/float64(ms_per_second))
	fmt.Println("stim_offset ", stim_offset/float64(ms_per_second))
	fmt.Println("stim_decay ", stim_decay/float64(ms_per_second))
	fmt.Println("odor_slip ", odor_slip)
	fmt.Println("odor_PN_set ", odor_PN_set)
	fmt.Println("odor_LN_set ", odor_LN_set)
	fmt.Println("BG_input_rate ", BG_input_rate)
	fmt.Println("ORN_input_rate ", ORN_input_rate)
	fmt.Println("BG_input_strength_PN ", BG_input_strength_PN)
	fmt.Println("BG_input_strength_LN ", BG_input_strength_LN)
	fmt.Println("ORN_input_strength_PN ", ORN_input_strength_PN)
	fmt.Println("ORN_input_strength_LN ", ORN_input_strength_LN)
	fmt.Println()
	// ...
	fmt.Println("time_end ", time_end/ms_per_second)
	fmt.Println("click_per_ms ", click_per_ms)
	fmt.Println("V_rest_PN ", V_rest_PN)
	fmt.Println("V_rest_LN ", V_rest_LN)
	fmt.Println("LN2PN_slow_prob ", LN2PN_slow_prob)
	fmt.Println("LN2PN_GABA_prob ", LN2PN_GABA_prob)
	fmt.Println("LN2LN_GABA_prob ", LN2LN_GABA_prob)
	fmt.Println("PN2PN_nACH_prob ", PN2PN_nACH_prob)
	fmt.Println("PN2LN_nACH_prob ", PN2LN_nACH_prob)
	// need to check the reality??
	fmt.Println("if_readin_matrix ", if_readin_matrix)
	fmt.Println("if_slowGABA_overlap ", if_slowGABA_overlap)
	fmt.Println("fading_rate ", fading_rate)
	fmt.Println("ms_per_frame ", ms_per_frame)
	fmt.Println()
	// ...
	fmt.Println("use_exist_config ", use_exist_config)
	fmt.Println("plot_file_dir ", plot_file_dir)
	fmt.Println("save_file_dir ", save_file_dir)
	fmt.Println("if_running ", if_running)
	fmt.Println("if_random_stim ", if_random_stim)
	fmt.Println("if_init_rt_plots ", if_init_rt_plots)
	fmt.Println("if_realtime_plot ", if_realtime_plot)
	fmt.Println("if_PN_plot ", if_PN_plot)
	fmt.Println("if_LN_plot ", if_LN_plot)
	fmt.Println("if_stimRaster_plot ", if_stimRaster_plot)
	fmt.Println("if_vRaster_plot ", if_vRaster_plot)
	fmt.Println("if_synaRaster_plot ", if_synaRaster_plot)
	fmt.Println("if_vFluct_plot ", if_vFluct_plot)
	fmt.Println("if_stimFluct_plot ", if_stimFluct_plot)
	fmt.Println("if_synaFluct_plot ", if_synaFluct_plot)
	fmt.Println()
}

func disp_spike_freq() {
	var i int64
	defer Exit(Enter("$FN()"))
	PN_spike_freq_file, err = os.OpenFile(PN_spike_freq_filename, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	LN_spike_freq_file, err = os.OpenFile(LN_spike_freq_filename, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
	check(err)
	// ...
	fmt.Println("\nSpiking rates (spikes per second) of each neuron <<", freq_scoping)
	for i = 0; i < PN_number; i++ {
		fmt.Printf("\t\tPN_spike_freq[%03d] \t %.1f\n", i, PN_spike_freq[i])
		_, _ = PN_spike_freq_file.Write([]byte(strconv.FormatInt(i, 10)))
		_, _ = PN_spike_freq_file.WriteString("   ")
		_, _ = PN_spike_freq_file.Write([]byte(strconv.FormatFloat(PN_spike_freq[i], 'f', 6, 64)))
		_, _ = PN_spike_freq_file.WriteString("\n")
	}
	// ...
	fmt.Println()
	for i = 0; i < LN_number; i++ {
		fmt.Printf("\t\tLN_spike_freq[%03d] \t %.1f\n", i, LN_spike_freq[i])
		_, _ = LN_spike_freq_file.Write([]byte(strconv.FormatInt(i, 10)))
		_, _ = LN_spike_freq_file.WriteString("   ")
		_, _ = LN_spike_freq_file.Write([]byte(strconv.FormatFloat(LN_spike_freq[i], 'f', 6, 64)))
		_, _ = LN_spike_freq_file.WriteString("\n")
	}
	// ...
	PN_spike_freq_file.Close()
	LN_spike_freq_file.Close()
	fmt.Println("Spiking rates saved to files.")
}

func init_odor() {
	var i int64
	defer Exit(Enter("$FN()"))
	odor_PN_set = mapset.NewSet()
	odor_LN_set = mapset.NewSet()
	for i = 0; i < PN_number; i++ {
		if i >= odor_slip && i < stim_PN_num+odor_slip {
			odor_PN_set.Add(i)
		}
	}
	for i = 0; i < stim_LN_num; i++ {
		odor_LN_set.Add(i)
	}
}

func proc_args() {
	defer Exit(Enter("$FN()"))
	arg_num := int64(len(os.Args))
	if arg_num > 3 {
		usage()
		return
	} else if arg_num == 3 { // reading user's odor_slip and matrix_seed
		t_int1, err1 := strconv.Atoi(os.Args[1])
		t_int2, err2 := strconv.Atoi(os.Args[2])
		if err1 != nil || err2 != nil { // error proc
			usage()
			check(err)
			return
		}
		// ...
		if t_int1 < 0 || int64(t_int1) >= PN_number {
			odor_slip = odor_slip % PN_number
			fmt.Println("WARNing: odor_slip is reset to ", odor_slip)
		} else {
			odor_slip = int64(t_int1)
		}
		matrix_seed = int64(t_int2)
		fmt.Println("using usr's odor_slip and matrix_seed")
	} else if arg_num == 2 { // reading user's odor_slip
		t_int, err := strconv.Atoi(os.Args[1])
		if err != nil { // error proc
			usage()
			check(err)
			return
		}
		// ...
		if t_int < 0 || int64(t_int) >= PN_number {
			odor_slip = odor_slip % PN_number
			fmt.Println("WARNing: odor_slip is reset to ", odor_slip)
		} else {
			odor_slip = int64(t_int)
		}
		fmt.Println("using usr's odor_slip and default(0) matrix_seed")
	} else if arg_num == 1 {
		fmt.Println("using default(0) odor_slip and matrix_seed")
	}
	// ...
	fmt.Println(odor_slip, " <= stim_PN_ID < ", stim_PN_num+odor_slip)
	fmt.Println("initial matrix with random seed: ", matrix_seed)
	fmt.Println()
}

func main() {
	var click int64
	var yaml_obj map[interface{}]interface{}
	var flagCpuprofile string
	var cpuProf_file *os.File
	//...
	defer Exit(Enter("$FN()"))
	if if_runtime_trace > 0 {
		flagCpuprofile = "cpuRuntime.prof"
		cpuProf_file, err = os.OpenFile(flagCpuprofile, os.O_RDWR|os.O_CREATE|os.O_TRUNC, 0640)
		check(err)
		pprof.StartCPUProfile(cpuProf_file)
		defer pprof.StopCPUProfile()
	}
	//...
	fmt.Println(time.Unix(time.Now().Unix(), 0).String())
	fmt.Println("Init...")
	proc_args()
	init_para()
	yaml_obj = decode_yaml(yaml_obj)
	proc_para(yaml_obj)
	init_matrix(yaml_obj)
	disp_matrix()
	save_matrix()
	init_rand()
	init_odor()
	//...
	fmt.Println("--- scale")
	scale_synapse_currents()
	remove_coupling_matrix()
	disp_para()
	//...
	fmt.Println("--- setup")
	init_neuron()
	proc_filenames()
	if if_init_rt_plots > 0 {
		init_plot()
	}
	//...
	fmt.Println("\nExec...")
	for click = 1; click <= click_num_total; click++ { // each time step is 0.01 ms
		//...
		if (click-1)%10000 == 9999 { // 0.1s
			yaml_obj = decode_yaml(yaml_obj) // reread in the configuration parameters
			proc_para(yaml_obj)
			fmt.Print("  ", float64(click)/100000.0)
		}
		//...
		CLOCK = CLOCK + time_stepLen
		take_iteration(click) // move a step forward
		//...
		if (click-1)%1000 == 999 && if_init_rt_plots > 0 && if_realtime_plot > 0 { // replot each 10 ms
			proc_plot_files(click / 100) // openFiles, write, closeFiles
			do_plot(click / 100)         // this is the main function for plotings
		}
		//...
		if if_running <= 0 { // paused by user
			yaml_obj = decode_yaml(yaml_obj)
			proc_para(yaml_obj) // read in the configuration parameters
			time.Sleep(10 * time.Millisecond)
		}
	}
	//...
	fmt.Println("\nDest...")
	if if_init_rt_plots > 0 {
		dest_plot()
	}
	disp_spike_freq() // disp_spike_freq() should be in ahead of disp_presync_num()!!!
	save_docs()
	//...
	fmt.Println("\nDone...")
	fmt.Println(time.Unix(time.Now().Unix(), 0).String())
}
