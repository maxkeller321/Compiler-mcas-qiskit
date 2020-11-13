import numpy as np 
import itertools
import os 
import sys 

sys.path.append("D:\Python\pi3diamond\pym8190a")
import pym8190a

from pi3diamond import pi3d
import UserScripts.helpers.snippets_awg as sna; reload(sna)
import AWG_M8190A_Elements as E


nuclei_list = ['14n', '13c414', '13c90']

transition_list = ['14n+1 ms0', '14n-1 ms0', '14n+1 ms-1', '14n-1 ms-1', '14n+1 ms+1', '14n-1 ms+1',
                   '13c414 ms0', '13c414 ms-1', '13c414 ms+1',
                   '13c90 ms0', '13c90 ms-1', '13c90 ms+1']

state_list = ['+++', '++-', '+-+', '+--', '0++', '0+-',
                    '0-+', '0--', '-++', '-+-', '--+', '---',
                    'n+', 'n-', 'nn+', 'nn-', '+', '-', '0']

single_nuclei_state_list = ['n+', 'nn+', '+', '-', '0']

def rz_nitrogen(mcas, theta, ms=0, mn=1, amp=1.0): 
    """
        rz rotation applied on the nitrogen nuclear spin. The sublevels mn = 0, 1 should be used for computation & the electron sublevels ms = 0, -1 
        params: 
            mcas: instance of the Multi-channel-sequence class
            theta: rotation angle in rad must be in the range of -2pi <= theta <= 2pi
            ms: electron spin state 
            mn: nuclear spins state: can just be mn =+/-1: mn =+1 stands for the transition between mn = 0 & mn = +1 and  mn = -1 for the transition between mn = 0 & mn = -1 
            amp: value between [0, 1] (defines the gate duration) A high amplitude of 1 leads to the fastest gate execution
    """
    ry_nitrogen(mcas, np.pi/2, ms, mn, amp)
    rx_nitrogen(mcas, theta, ms, mn, amp)
    ry_nitrogen(mcas, -np.pi/2, ms, mn, amp)

def rx_nitrogen(mcas, theta, ms=0, mn=1, amp=1.0):
    """
        rx rotation applied on the nitrogen nuclear spin. The sublevels mn = 0, 1 should be used for computation  & the electron sublevels ms = 0, -1 
        params: 
            mcas: instance of the Multi-channel-sequence class
            theta: rotation angle in rad must be in the range of -2pi <= theta <= 2pi
            ms: electron spin state 
            mn: nuclear spins state: can just be mn =+/-1: mn =+1 stands for the transition between mn = 0 & mn = +1 and  mn = -1 for the transition between mn = 0 & mn = -1 
            amp: value between [0, 1] (defines the gate duration) A high amplitude of 1 leads to the fastest gate execution
    """
    transition = get_transition('14n', ms=ms, mn=mn)
    nuclear_rotation(mcas, theta, 'x', transition, amp)

def ry_nitrogen(mcas, theta, ms=0, mn=1, amp=1.0):
    """
        ry rotation applied on the nitrogen nuclear spin.The sublevels mn = 0, 1 should be used for computation  & the electron sublevels ms = 0, -1 
        params: 
            mcas: instance of the Multi-channel-sequence class
            theta: rotation angle in rad must be in the range of -2pi <= theta <= 2pi
            ms: electron spin state 
            mn: nuclear spins state: can just be mn =+/-1: mn =+1 stands for the transition between mn = 0 & mn = +1 and  mn = -1 for the transition between mn = 0 & mn = -1 
            amp: value between [0, 1] (defines the gate duration) A high amplitude of 1 leads to the fastest gate execution
    """
    transition = get_transition('14n', ms=ms, mn=mn)
    nuclear_rotation(mcas, theta, 'y', transition, amp)


def rz_carbon_90(mcas, theta, ms=-1, amp=1.0):
    """
        rz rotation applied on the carbon 90 nuclear spin. The electron sublevels ms = 0, -1 are used for computation 
        params: 
            mcas: instance of the Multi-channel-sequence class
            theta: rotation angle in rad must be in the range of -2pi <= theta <= 2pi
            ms: electron spin state 
            amp: value between [0, 1] (defines the gate duration) A high amplitude of 1 leads to the fastest gate execution
    """
    ry_carbon_90(mcas, np.pi/2, ms, amp)
    rx_carbon_90(mcas, theta, ms, amp)
    ry_carbon_90(mcas, -np.pi/2, ms, amp)

def rx_carbon_90(mcas, theta, ms=-1, amp=1.0):
    """
        rx rotation applied on the carbon 90 nuclear spin. The electron sublevels ms = 0, -1 are used for computation 
        params: 
            mcas: instance of the Multi-channel-sequence class
            theta: rotation angle in rad must be in the range of -2pi <= theta <= 2pi
            ms: electron spin state 
            amp: value between [0, 1] (defines the gate duration) A high amplitude of 1 leads to the fastest gate execution
    """
    transition = get_transition('13c90', ms=ms)
    nuclear_rotation(mcas, theta, 'x', transition, amp)

def ry_carbon_90(mcas, theta, ms=-1, amp=1):
    """
        ry rotation applied on the carbon 90 nuclear spin.  The electron sublevels ms = 0, -1 are used for computation 
        params: 
            mcas: instance of the Multi-channel-sequence class
            theta: rotation angle in rad must be in the range of -2pi <= theta <= 2pi
            ms: electron spin state 
            amp: value between [0, 1] (defines the gate duration) A high amplitude of 1 leads to the fastest gate execution
    """
    transition = get_transition('13c90', ms=ms)
    nuclear_rotation(mcas, theta, 'y', transition, amp)

def rz_carbon_414(mcas, theta, ms=-1, amp=1.0):
    """
        rz rotation applied on the carbon 414 nuclear spin. The electron sublevels ms = 0, -1 are used for computation 
        params: 
            mcas: instance of the Multi-channel-sequence class
            
            theta: rotation angle in rad must be in the range of -2pi <= theta <= 2pi
            ms: electron spin state 
            amp: value between [0, 1] (defines the gate duration) A high amplitude of 1 leads to the fastest gate execution
    """
    ry_carbon_414(mcas, np.pi/2, ms, amp)
    rx_carbon_414(mcas, theta, ms, amp)
    ry_carbon_414(mcas, -np.pi/2, ms, amp)

def rx_carbon_414(mcas, theta, ms=-1, amp=1.0):
    """
        rx rotation applied on the carbon 414 nuclear spin.  The electron sublevels ms = 0, -1 are used for computation 
        params: 
            mcas: instance of the Multi-channel-sequence class
            theta: rotation angle in rad must be in the range of -2pi <= theta <= 2pi
            ms: electron spin state 
            amp: value between [0, 1] (defines the gate duration) A high amplitude of 1 leads to the fastest gate execution
    """
    transition = get_transition('13c414', ms=ms)
    nuclear_rotation(mcas, theta, 'x', transition, amp)

def ry_carbon_414(mcas, theta, ms=-1, amp=1.0):
    """
        ry rotation applied on the carbon 414 nuclear spin.  The electron sublevels ms = 0, -1 are used for computation 
        params: 
            mcas: instance of the Multi-channel-sequence class
            theta: rotation angle in rad must be in the range of -2pi <= theta <= 2pi
            ms: electron spin state 
            amp: value between [0, 1] (defines the gate duration) A high amplitude of 1 leads to the fastest gate execution
    """
    transition = get_transition('13c414', ms=ms)
    nuclear_rotation(mcas, theta, 'y', transition, amp)

def rz_electron(mcas, theta, nuclear_spin_state, amp=1.0, mixer_deg=-90):
    electron_rotation(mcas,  np.pi/2, 'y', nuclear_spin_state, amp=amp, mixer_deg=mixer_deg)
    electron_rotation(mcas, theta, 'x', nuclear_spin_state, amp=amp, mixer_deg=mixer_deg)
    electron_rotation(mcas, -np.pi/2, 'y', nuclear_spin_state, amp=amp, mixer_deg=mixer_deg)

def rx_electron(mcas, theta, nuclear_spin_state, amp=1.0, mixer_deg=-90):
    electron_rotation(mcas, theta, 'x', nuclear_spin_state, amp=amp, mixer_deg=mixer_deg)

def ry_electron(mcas, theta, nuclear_spin_state, amp=1.0, mixer_deg=-90):
    electron_rotation(mcas, theta, 'y', nuclear_spin_state, amp=amp, mixer_deg=mixer_deg)

def nuclear_rotation(mcas, theta, rotation_axis, transition, amp):
    """
        rx rotation applied on the nitrogen nuclear spin. The sublevels mn = 0, 1 are used for computation & the electron sublevels ms = 0, -1 
        params: 
            mcas: instance of the Multi-channel-sequence class
            theta: rotation angle in rad must be in the range  of -2pi <= theta <= 2pi
            rotation_axis: axis of rotation 'x' or 'y' (to make an rotation around 'z' use another function)
            transition: key from get_transition that describes the correct rotation 
            amp: value between [0, 1] (defines the gate duration) A high amplitude of 1 leads to the fastest gates execution
    """

    if transition not in transition_list: 
        raise Exception("The Transition is not a valid transition")

    theta, phase = get_optimised_angle_and_phase(theta, rotation_axis)

    scaling_factor = theta/np.pi 
    lenth_mus =  E.round_length_mus_full_sample(pi3d.tt.rp(transition, amp=amp).pi*scaling_factor) 

    sna.nuclear_rabi(mcas,
                    name=transition,
                    frequencies=[pi3d.tt.t(transition).current_frequency],
                    amplitudes=[amp],
                    phases=[phase],
                    length_mus=lenth_mus,
                    new_segment=True)
    mcas.asc(length_mus=10) # delay for RF field fluctuations

def electron_rotation(mcas, theta, rotation_axis, nuclear_spin_state, amp=1.0, mixer_deg=-90):
    """
        This function is not yet tested! 

        This function applies a non robust electron rotation around the angle theta for a certain nuclear spin state
        params: 
            mcas: instance of the Multi-channel-sequence class
            theta: rotation angle in rad must be in the range  of -2pi <= theta <= 2pi
            rotation_axis: axis of rotation 'x' or 'y' (to make an rotation around 'z' use another function)
            nuclear_spin_state: certain nuclear spin stat must be one of: ['+++', '++-', '+-+', '+--', '0++', '0+-',
                    '0-+', '0--', '-++', '-+-', '--+', '---']
            amp: value between [0, 1] (defines the gate duration) A high amplitude of 1 leads to the fastest gates execution
            mixer_deg: 
    """
    # at the moment this functions is just build for a fully known register and not for +,-,0 and n+ and nn+
    theta, phase = get_optimised_angle_and_phase(theta, rotation_axis)

    freq = pi3d.tt.mfl({'14N': [{'+': +1, '0': 0, '-': -1}[nuclear_spin_state[0]]], '13c414': [{'+': +.5, '-': -.5}[nuclear_spin_state[1]]], '13c90': [{'+': +.5, '-': -.5}[nuclear_spin_state[2]]]}, ms_trans='left')
    scaling_factor = theta/np.pi 
    lenth_mus=pi3d.tt.rp('e_rabi', mixer_deg=mixer_deg, amp=amp).pi*scaling_factor
    lenth_mus = E.round_length_mus_full_sample(lenth_mus)
    sna.electron_rabi(mcas,
                        name='electron rabi',
                        length_mus=lenth_mus,
                        amplitudes=[amp],
                        frequencies=freq,
                        phases=[phase],
                        new_segment=True,
                        mixer_deg=mixer_deg)

def get_optimised_angle_and_phase(theta, rotation_axis):
    """
        optimises the angle theta (e.g rotation theta > np.pi --> makes more sense to go for a negative rotation )
        and returns the correct phase for the AWG
        params: 
            theta: rotation angle in rad must be in the range of -2pi <= theta <= 2pi
            rotation_axis: axis of rotation 'x' or 'y' (to make an rotation around 'z' use another function)
    """

    if (type(theta) == int) or (type(theta) == float):   
        if (-2*np.pi <= theta < -np.pi): 
            postitive_rotation = True 
            theta += 2*np.pi
        elif (-np.pi <= theta < 0):
            postitive_rotation = False
            theta = abs(theta) 
        elif (0 <= theta <= np.pi):
            postitive_rotation = True 
            theta = theta
        elif (np.pi < theta <= 2*np.pi): 
            postitive_rotation = False
            theta = 2*np.pi - theta
        else: 
            raise Exception("The rotation angel is not in the range of -2pi <= theta <= 2pi. Insert a valid radian value!")
    else: 
        raise Exception("The rotation angle theta must be from the type integer or float")


    if str(rotation_axis).strip() == 'x': 
        if postitive_rotation: 
            phase = 0.0 # phase for the waveform generator in degree 
        else: 
            phase = 180.0 

    elif str(rotation_axis).strip() == 'y':
        if postitive_rotation:
            phase = 90.0  
        else: 
            phase = -90.0  
    else: 
        raise Exception("Only 'x' and 'y' are valid rotation axis!")

    return theta, phase


def get_transition(nucleus, ms=0, mn=1):
    """
        returns the correct transition for the rabis based on the passed parameters: 
        One from: ['14n+1 ms0', '14n-1 ms0', '14n+1 ms-1', '14n-1 ms-1', '14n+1 ms+1', '14n-1 ms+1',
                '13c414 ms0', '13c414 ms-1', '13c414 ms+1',
                '13c90 ms-1', '13c90 ms+1']
        params: 
            nucleus: "14n", '13c414' or '13c90'
            ms: electron spin state in the moment: -1, 0, +1
            mn: nitrogen nuclear spins state: can just be mn =+/-1: mn =+1 stands for the transition between mn = 0 & mn = +1 and 
             mn = -1 for the transition between mn = 0 & mn = -1 
    """

    if "14n" in nucleus.lower():
        if mn == -1:
            nuc_part = '14n-1' 
        elif mn == 1: 
            nuc_part = '14n+1' 
        else: 
            raise Exception("Only mn = -1 or mn = 1 is allowed!")
    elif '13c414' in nucleus.lower():
        nuc_part = '13c414' 
    elif '13c90' in nucleus.lower(): 
        nuc_part = '13c90' 
    else: 
        raise Exception("Inserted nucleus is not valid!")
    
    if ms == 1: 
        electron_part = 'ms+1'
    elif ms == 0: 
        electron_part = 'ms0'
    elif ms == -1: 
        electron_part = 'ms-1'
    else: 
        raise Exception("Only ms = +1, 0 or ms = -1 is allowed!")

    transition = '{} {}'.format(nuc_part, electron_part) 

    if transition not in transition_list: 
        raise Exception("The Transition: {} is not a valid transition".format(transition))

    return transition


def readout_nuclear_spin_state(mcas, state, step_idx=0): 
    """
        Appends a single shot readout sequence of a certain nuclear spin state to the passed mcas.
        params: 
            mcas: instance of the Multi-channel-sequence class
            state: state which should be read out (which controlls the electron rabi)
                Pick one from:  ['+++', '++-', '+-+', '+--', '0++', '0+-',
                            '0-+', '0--', '-++', '-+-', '--+', '---',
                                'n+', 'nn+', '+', '-', '0']
            step_idx: If multiple SSR's are used, the step_idx must refer to the corresponding ssr 
                e.g.: (first ssr: step_idx = 0, second ssr: step_idx = 1 ...)
    """

    def get_nuclei(state):

        if type(state) == str:
            state = state.strip().lower()
            if state =='+': 
                return '14n+1' 
            elif state == '0': 
                return '14n0' 
            elif state == '-': 
                return '14n-1' 
            elif state == 'n+' or state == 'n-':
                return '13c414'
            elif state == 'nn+' or state == 'nn-':
                return '13c90'
            elif state in ["".join(i) for i in itertools.product(['+', '0', '-'], ['+', '-'], ['+', '-'])]: 
                return '13c90'
            else: 
                raise Exception("The entered state is not valid. State must be one of: \n {}".format(state_list))  

        else: 
            raise Exception("Invalid state datatype! state must be of type string!")
    
    def get_inverse_state(state): 
        if '+' in state and state.count('n') > 0: 
            return state.replace('+', '-')
        elif '-' in state and state.count('n') > 0: 
            return state.replace('-', '+')
        else: 
            raise Exception("Inserted state is not correct")
        
    nucleus = get_nuclei(state)
    if (nucleus == '13c90' and state.count('n') == 2) or nucleus == '13c414':

        wave_file_kwargs_all_but_standard = dict(filepath=sna.wfpd_standard[get_inverse_state(state)], rp=pi3d.tt.rp('e_rabi', mixer_deg=-90))
        wave_file_kwargs_all = dict(filepath=sna.wfpd_standard[state], rp=pi3d.tt.rp('e_rabi', mixer_deg=-90))

        sna.ssr(mcas,
                transition='left',
                robust=True, # necessary
                laser_dur=sna.__LASER_DUR_DICT__[nucleus],
                mixer_deg=-90,
                nuc=nucleus,
                frequencies=[pi3d.tt.mfl({'14n': [0]}), pi3d.tt.mfl({'14n': [0]})], # frequencies are not used, when we pass a wavefile to the ssr
                wave_file_kwargs=[wave_file_kwargs_all, wave_file_kwargs_all_but_standard],
                repetitions=sna.__SSR_REPETITIONS__[nucleus],
                step_idx=step_idx)

    else:
        sna.ssr_single_state(mcas, state, repetitions=sna.__SSR_REPETITIONS__[nucleus], step_idx=step_idx)
        

def initialise_nuclear_spin_register(mcas, state):
    """
        Appends an initialisation sequence to the passed mcas.
        params: 
            mcas: instance of the Multi-channel-sequence class
            state: state in which the register should be initialised 
                Pick one from:  ['+++', '++-', '+-+', '+--', '0++', '0+-',
                            '0-+', '0--', '-++', '-+-', '--+', '---',
                                'n+', 'nn+', '+', '-', '0']
    """
    if type(state) == str:
        state = state.strip().lower()
        if state in state_list:
            if state in ['+', '-', '0']:
                sna.init_14n(mcas, new_segment=True, mn=state)
            elif state in ['n+', 'n-']:
                sna.init_13c(mcas, s='414', new_segment=True, state=state[-1])
            elif state in ['nn+', 'nn-']: 
                sna.init_13c(mcas, s='90', new_segment=True, state=state[-1])
            else: 
                sna.init_14n(mcas, new_segment=True, mn=state[0])
                sna.init_13c(mcas, s='414', new_segment=False, state=state[1])
                sna.init_13c(mcas, s='90', new_segment=False, state=state[2])
        else:
            raise Exception("The entered state is not valid. State must be one of: \n {}".format(state_list))  
    else:
        raise Exception("Invalid state datatype! state must be of type string!")

def initialise_electron_spin(mcas): 
    """
        Append an electron spin initialisation into ms=0
        Params: 
            mcas: instance of the Multi-channel-sequence class
    """
    sna.polarize_green(mcas, new_segment=True)

def initialise_with_red(mcas):
    """
            Append an electron spin initialisation into ms=0
        Params:
            mcas: instance of the Multi-channel-sequence class
    """
    sna.polarize(mcas, new_segment=True)


def cnot_between_nuclear_spins(mcas, controlled_qubit, controlling_qubit): 
    """
        Not yet tested! This function should append a controlled not gate to the mcas between two nuclei!  
        params: 
            mcas: instance of the Multi-channel-sequence class
            controlled_qubit: controlled nuclear spin: (14n, 13c414, or 13c90) 
            controlling_qubit: controlling nuclear spin: (14n, 13c414, or 13c90) 
    """

    if type(controlling_qubit) != str or type(controlled_qubit) != str: 
        raise Exception("Controlled & controlling nuclei must be from the type string!") 
    if controlled_qubit.lower() not in nuclei_list: 
        raise Exception("Controlled nuclei not valid! Must be one of: \n {}".format(nuclei_list))
    if  controlling_qubit.lower() not in nuclei_list: 
        raise Exception("Controlling nuclei not valid! Must be one of: \n {}".format(nuclei_list))

    initialise_electron_spin(mcas)

    if controlling_qubit.lower() == '14n': 
        electron_controlled_not(mcas, '+')
    elif controlling_qubit.lower() == '13c414': 
        electron_controlled_not(mcas, 'n+')
    elif controlling_qubit.lower() == '13c90': 
        electron_controlled_not(mcas, 'nn+')

    if controlled_qubit.lower() == '14n': 
        rx_nitrogen(mcas, np.pi, ms=-1)
    elif controlled_qubit.lower() == '13c414': 
        rx_carbon_414(mcas, np.pi, ms=-1)
    elif controlled_qubit.lower() == '13c90':
        rx_carbon_90(mcas, np.pi, ms=-1)

def electron_pi_pulse(mcas, ms_transition='left'):
    """
        not fully written 
        & not tested 
        params: 
            mcas: instance of the Multi-channel-sequence class
            ms_transition: stands for the addressed electron sub levels 
                            ('left' stands for transition between ms = -1/ 0 and 'right' for the transition between ms=+1/0 )

    """
    
    sna.single_robust_electron_pi(mcas,
                                    nuc='all',
                                    transition=ms_transition, # could also be -1 (needs to be tested)
                                    frequencies=pi3d.tt.mfl({'14n': [0]}, ms_trans={'left':'-1', 'right':'+1'}[ms_transition]),
                                    new_segment=True)


def electron_controlled_not(mcas, state): 
    """
        This function is not yet tested. There are two implementations in the function one should be commented out!
            Both need to be tested. If the implementation with the electron_rabis is working then this implementation 
            should be prefered!
        
        It uses the optimal control pulses from snippets to realize 
        controlled not gates on the electron spin dependent on a certain nuclear spin state. 
        params: 
            mcas: instance of the Multi-channel-sequence class
            state: nuclear spin state on which the electron rotation should be dependend 
                Pick one from:  ['+++', '++-', '+-+', '+--', '0++', '0+-',
                            '0-+', '0--', '-++', '-+-', '--+', '---',
                                'n+', 'nn+', '+', '-', '0']
    """
    """
    if state not in state_list: 
        raise Exception("The entered state is not valid. State must be one of: \n {}".format(state_list))  

    
    mcas.start_new_segment(name="cx", loop_count=1)

    wave_file_kwargs = [dict(
                    filepath=sna.wfpd_standard[state],
                    rp=pi3d.tt.rp('e_rabi', mixer_deg=-90)
                    )]
    wave_file = [E.WaveFile(**i) for i in wave_file_kwargs][0]
    length_mus_mw = wave_file.length_mus


    def part_step():
        t = np.cumsum([0.0, length_mus_mw])
        return {
            2: wave_file.ret_part(t[0], t[1]),
        }

    def pd2g_dict(mixer_deg):

        pd2g_dict = {1: {}, 2: {}}
        for ch in [1, 2]:
            # for i in [2, 3, 4]:
            for i in [2]:
                pd2g_dict[ch][i] = {}
                pd2g_dict[ch][i]['frequencies'] = pi3d.tt.mfl({'14n': [0]}) # is not going to be used but needs to be passed
                pd2g_dict[ch][i].update(dict(type='robust', wave_file=E.WaveFile(part=part_step()[i], **wave_file_kwargs[0])))
                if ch == 2:
                    pd2g_dict[ch][i]['phases'] = np.array([mixer_deg])
                else:
                    pd2g_dict[2][i] = {}

        return pd2g_dict
    
    d = pd2g_dict(-90.0)
    aa = dict() # I don't think this is needed but I'm still going to pass it
    mcas.asc(pd2g1=d[1][2], pd2g2=d[2][2], name='MW', **aa)
    """

    """Second implementation"""


    wave_file = E.WaveFile(filepath=sna.wfpd_standard[state],
                            nonlinear_params=pi3d.tt.e_rabi_left.process_electron_rabi_file()) # maybe consider also right


    sna.electron_rabi(mcas,
                        #name='electron rabi',
                        frequencies=pi3d.tt.mfl({'14N': [0]}),
                        wave_file=wave_file,
                        wait_switch=0.0,
                        new_segment=True)


def electron_controlled_not_2(mcas, nucstate):

    """
    Not yet tested & it should only be used when electron_controlled_not is not working!  
    Controlled not gate on the electron spin dependent on the nuclear spin state
    There is a documentation for nucstate in one of the oneNote pages 
    """
    cnot_dict={}
    for root, dirs, files in os.walk(r"\\PI3-PC161\d\Python\pi3diamond\UserScripts\Robust\new_pulses_for_qutrit_qft_sensing", topdown=False):
        for name in files:
            if name in ['mw_aphi.dat', 'MW.dat']:
                cnot_dict[os.path.basename(root)] = os.path.join(root, name)

    for key, val in cnot_dict.items():
        if 'FN' in key:
            k = key[::-1][key[::-1].index('NF')+2:]
            cnot_dict[k[:k.index('_')][::-1]] = val
            del cnot_dict[key]
        if '_qutrit_cphase_' in key:
            k = key[len('20170317-h18m07s08_qutrit_'):]
            cnot_dict[k] = val
            del cnot_dict[key]

    wfd = {}
    for key, val in cnot_dict.items():
        wfd[key] = pym8190a.elements.WaveFile(filepath=val, rp=pi3d.tt.rp('e_rabi', mixer_deg=-90))



    def crot(nucstate):
        sna.electron_rabi(mcas,
                            name='electron pi',
                            wave_file=wfd[nucstate],
                            frequencies=pi3d.tt.mfl({'14N': [0]}),
                            new_segment=False,
                            mixer_deg=-90)

    crot(nucstate)