""" Copyright © 2023, Spencer Hiscox. All rights Reserved.
"""


from time import perf_counter as pctime
from time import sleep
from math import sqrt, pi, sin, cos, tan
π = pi
from typing import Tuple
import re
from gc import collect
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from copy import deepcopy



#CLASS DEFINITIONS

class truss:
    def __init__(self, no_joists, angle, t_type):
        self.no_joists = no_joists
        self.m_dimensions = {'angle': angle, 'Sj': _span / (no_joists + 1)}
        self.m_dimensions['Diag'] = self.m_dimensions['Sj'] / angle[1]
        self.m_dimensions['Vert'] = self.m_dimensions['Diag'] * angle[2]
        self.m_weights = [[], [], [], [], []]
        self.m_forces = [[], [], [], [], []]
        self.m_sections = [[], [], [], [], []]
        self.total_weight = 0
        
        self.typ = t_type
        
        if self.typ == 3:
            self.m_dimensions['Horz'] = self.m_dimensions['Sj'] * 2
            self.m_dimensions['Bot'] = self.m_dimensions['Horz'] * 2
        if self.typ == 2:
            self.m_dimensions['Bot'] = self.m_dimensions['Sj'] * 2



class diag_obj:
    
    def __init__(self):
        pass
    
    def forces(self, **kwargs)->None:
        
        global t_count, t0
        n_upd = kwargs.get('no_update')
        user_choice = []
        user_choice.append(input(f"{new_screen}\nView graphs from the two algorithms' generated data [Y / N]? ").lower())
        if user_choice[0] == 'y' or user_choice == 'yes':
            user_choice[0] == True
        else:
            user_choice[0] = False
        user_choice.append(input("Display Optimal Truss solution upon completion of optimization algorithms [Y / N]?  ").lower())
        if user_choice[1] == 'y' or user_choice == 'yes':
            user_choice[1] = True
        else:
            user_choice[1] = False
        t_count = 0
        t0 = pctime()
        optimize(disp_res=user_choice[1], disp_graphs=user_choice[0])
        retainer = [weights_by_type, deepcopy(angle_matrix)]
        input("Fast Optimization complete. Press <ENTER> to Continue...")
        t_count = 0
        t0 = pctime()
        optimize(basic=True, disp_res=user_choice[1], disp_graphs=user_choice[0])
        input("Hard-code calculation of forces complete. Press <ENTER> to compare algorithm results...")
        print(f"{new_screen}\nVerifying all weights match, comparing by angle and truss type...\n\n")
        sleep(3)
        fail = False
        truss_types = [f"{'Howe':>24}\t", f"{'Pratt':>24}\t", f"{'Warren with Verticals':>24}\t", f"{'Warren without Verticals':>24}\t"]
        for truss_type in range(len(weights_by_type)):
            row_no, row_change = 1, False
            print(truss_types[truss_type], end="")
            for angle_entry in range(len(weights_by_type[truss_type])):
                if weights_by_type[truss_type][angle_entry] != retainer[0][truss_type][angle_entry]:
                    fail = True
                if row_change:
                    print(f"{('False', 'True')[not fail]:>36}", end="\t")
                    row_change = False
                else:
                    print(not fail, end="\t")
                if angle_entry > (131 / 8) * row_no - 1:
                    row_no += 1
                    row_change = True
                    print()
                sleep(0.1)
            print()
        if not fail:
            if not n_upd:
                update_truss_matrix(retainer[1])
            print("\n\nALL CALCULATED VALUES MATCH. FAST ALGORITHM IS VERIFIED. DIAGNOSTICS COMPLETE.\n")
            input("Press <ENTER> to Continue...")
        else:
            input("\n\nTEST FAILED. FAST ALGORITHM IS NON-FUNCTIONAL. PLEASE USE HARD-CODE (BASIC) ALGORITHM FOR FORCE CALCULATIONS.\n"
                  "(HARD-CODE OPTIMIZATION RESULTS WILL BE STORED IN MEMORY WHEN YOU EXIT THESE DIAGNOSTICS)\n"
                  "DIAGNOSTICS COMPLETE.\nPress <ENTER> to Continue...")
            if not n_upd:
                update_truss_matrix()
        del retainer
        collect();        



#GLOBAL FUNCTION DEFINITIONS

def disp_truss_menu()->None:
    
    def invalid_response(string: str)->str:
        print(f"{string} is not a valid response.")
        input("Press <ENTER> to continue.")
        print(f"{new_screen}")
    
    if not truss_matrix[0]:
        input(f"{new_screen}No Optimization results loaded into memory. Run optimization or load saved results from file.\n"
              "Press <ENTER> to Continue...")
    else:
        user_choice = []
        while True:
            print(f"{new_screen}\t1:\tHowe\n\t2:\tPratt\n\t3:\tWarren (with verticals)\n\t4:\tWarren (without verticals)\n\n")
            user_choice.append(input("Select truss type from the menu options: ").lower())
            try: 
                user_choice[0] = int(user_choice[0])
                if (0 < user_choice[0] < 5):
                    user_choice[0] -= 1
                    break
                else:
                    invalid_response(user_choice[0])
                    continue
            except:
                try:
                    if re.findall('how', user_choice[0]):
                        user_choice[0] = 0
                    elif re.findall('prat', user_choice[0]):
                        user_choice[0] = 1
                    elif re.findall('warr?en', user_choice[0]):
                        if re.findall('Ø', user_choice[0]) or re.findall('no', user_choice[0]) or \
                           re.findall('nv', user_choice[0]):
                            user_choice[0] = 3
                        else:
                            user_choice[0] = 2
                    else:
                        invalid_repsonse(user_choice[0])
                        continue
                    break
                except:
                    invalid_response(user_choice[0])
                    continue
        
        while True:
            user_choice.append(input("\nEnter angle: ").lower())
            if user_choice[1].isdigit():
                user_choice[1] = int(user_choice[1])
                if user_choice [1] >= 90 or user_choice[1] <= 0:
                    invalid_response(user_choice[1])
                    continue
                break
            else:
                try:
                    user_choice[1] = float(user_choice[1])
                    if int(user_choice[1]) >= 90 or int(user_choice[1]) <= 0:
                        invalid_response(user_choice[1])
                        continue                    
                    break
                except:
                    invalid_response(user_choice[1])
                    continue
        try:
            display_truss(truss_matrix[user_choice[0]][user_choice[1]])
        except:
            angle = user_choice[1] * pi / 180
            angle = [angle, sin(angle), cos(angle), tan(angle)]
            truss_matrix[user_choice[0]][user_choice[1]] = (fast_calculate_member_forces(number_joists, angle, single=True, typ=user_choice[0]),
                                                                 calculate_member_forces(number_joists, angle, single=True, typ=user_choice[0]))[basic_force_algorithm]
            assign_HSS_sections(truss_matrix[user_choice[0]][user_choice[1]])
            display_truss(truss_matrix[user_choice[0]][user_choice[1]])



def display_truss(obj: '__main__.truss')->None:
    
    truss_names = ['Howe', 'Pratt', 'Warren\n\t\t(with vertical members)', 'Warren\n\t\t(without vertical members)']
    bullet = '● '
    horizontal_display = f'{obj.m_dimensions["Sj"]}\n'
    if obj.typ > 1:
        if obj.typ > 2:
            horizontal_display = f'\n{"GENERAL HORIZONTALS:":>36}\t{obj.m_dimensions["Horz"]}' + f'\n{"BOTTOM CENTER:":>36}\t{obj.m_dimensions["Bot"]}\n'
        else:
            horizontal_display = '\n\t\t{"TOP MEMBERS:":>36}\t' + horizontal_display + f'\n\t\t{"BOTTOM MEMBERS:":>36}\t{obj.m_dimensions["Bot"]}\n'
    row_headers = ['TOP MEMBERS', 'OUTWARD-LEANING DIAGONAL MEMBERS', 'VERTICAL MEMBERS', 'INWARD-LEANING DIAGONAL MEMBERS', 'BOTTOM MEMBERS']
    if obj.total_weight == None:
        print(f"{new_screen}TRUSS DESIGN / CONFIGURATION DISPLAY\n{'':=>36}\n")
        print("< UNABLE TO RESOLVE TRUSS DESIGN:  HSS SECTIONS AVAILABLE WILL FAIL UNDER OPERATIONAL LOADS >\n\n"
              f'{"TYPE:":>9}\t{truss_names[obj.typ]}\n'
              f'{"ANGLE*:":>9}\t{obj.m_dimensions["angle"][0] / pi * 180:.3f}\n'
              f'WEIGHT**:\t{obj.total_weight}\n\n\n'
              f'*ALL ANGLES ARE MEASURED TO THE NORMAL AND GIVEN IN DEGREES (°)\n**WEIGHT GIVEN IN KILONEWTONS (kN)\n{"":=>156}\n', sep="")        
    else:
        print(f"{new_screen}TRUSS DESIGN / CONFIGURATION DISPLAY\n{'':=>36}\n")
        print(((f"\t<OPTIMAL {truss_names[obj.typ]} TRUSS>\n\n", "\t<OPTIMAL DESIGN>\n\n")[obj == optimal_truss_design])[obj in optima] +
              f'{"TYPE:":>9}\t{truss_names[obj.typ]}\n'
              f'{"ANGLE*:":>9}\t{obj.m_dimensions["angle"][0] / pi * 180:.3f}\n'
              f'WEIGHT**:\t{obj.total_weight}\n\n\n'
              f'*ALL ANGLES ARE MEASURED TO THE NORMAL AND GIVEN IN DEGREES (°)\n**WEIGHT GIVEN IN KILONEWTONS (kN)\n{"":=>156}\n')
        print(f'DIMENSIONS*:\n{"":=>11}\n\n'
              f'\t  VERTICAL MEMBERS:\t{obj.m_dimensions["Vert"]}\n'
              f'\t  DIAGONAL MEMBERS:\t{obj.m_dimensions["Diag"]}\n'
              f'\tHORIZONTAL MEMBERS:\t{horizontal_display}\n\n'
              f'*ALL DIMENSIONS ARE GIVEN IN METERS (m)\n')
        print(f'NOTE:\n\t{bullet}THIS TRUSS DESIGN / CONFIGURATION IS PREDICATED ON A JOIST SPACING OF {obj.m_dimensions["Sj"]}m.\n'
              f'\t{bullet}THIS JOIST SPACING ALLOWS FOR {obj.no_joists} JOISTS EVENLY DISTRIBUTED ACROSS A ROOF SPAN OF {_span}m\n\n{"":=>156}\n\n')
        print(f'HSS SECTIONS, BY MEMBER*:\n{"":=>24}\n')
        for i in range(len(obj.m_sections)):
            print(f"{f'{row_headers[i]}':>32}\t", end="")
            for entry in obj.m_sections[i]:
                print(entry, end="\t")
            print()
        print(f"\n\n*HSS SECTIONS ARE LISTED BY MEMBER, IN ORDER FROM THE OUTERMOST MEMBER TO THE INNERMOST (CENTRAL) MEMBERS\nAS VIEWED, LEFT TO RIGHT IN THE ABOVE TABLE\n{'':=>156}\n\n")
        print(f"INTERNAL MEMBER FORCES*, BY MEMBER:\n{'':=>34}\n")
        for i in range(len(obj.m_forces)):
            print(f"{f'{row_headers[i]}':>32}", end="")
            for entry in obj.m_forces[i]:
                print(f"{f'{entry:.3f}':>16}", end="")
            print()
        print("\n\n*ALL FORCES ARE GIVEN IN KILONEWTONS (kN)")
    
    if disp_from_disp:
        while True:
            user_choice = input("\n[D]isplay another truss or Return to [M]ain Menu [D / M]? ").lower()
            if user_choice == 'd' or re.findall('display', user_choice):
                disp_truss_menu()
                break
            elif user_choice == 'm' or re.findall('main', user_choice) or re.findall('menu', user_choice):
                break
            else:
                input("Invalid response. Please enter either <D> or <M> to make your selection.\n"
                      "Press <ENTER> to Continue...")



def c_round(num: float, precision: int = 0)->int:
    
    if precision < 0:
        precision = abs(precision)
        subtraction = num % (10 ** precision)
        tmp = int(subtraction / 10 ** (precision - 1))
        return (int(num - subtraction), int(num - subtraction + 10 ** precision))[tmp > 4]
    else:
        #precision += 1
        tmp = int((num * 10 ** (precision + 1)) % 10)
        if precision == 0:
            return (int(num), int(num) + 1)[tmp > 4]
        else:
            return (int(num * 10 ** precision) / 10 ** precision, (int(num * 10 ** precision) + 1) / 10 ** precision)[tmp > 4]



def full_list(**kwargs)->None:
    
    from_update = kwargs.get('f_update')
    first_run = kwargs.get('first_rn')
    
    
    def left_fmt(*args)->str:
        if len(args) > 2:
            raise ValueError("Function left_fmt() takes a maximum of 2 arguments.")
        string, space = '', 75
        if len(args) == 2:
            if type(args[0]) == str and type(args[1]) == int:
                string, space = args[0], args[1]
            elif type(args[0]) == int and type(args[1]) == str:
                string, space = args[1], args[0]
            else:
                raise ValueError("Function left_fmt() can take only one argument of "
                                 "type string and one argument of type int.")
        if len(args) == 1:
            if type(args[0]) == str:
                string = args[0]
            elif type(args[0]) == int:
                space = args[0]
            else:
                raise ValueError("Function left_fmt() can take only arguments of "
                                 "type str or type int.")
        space -= len(string)
        space *= ' '
        return f"{string}{space}"
    
    headers = [f'{"F17 SPECIFIC VARIABLES:":<75}CONSTANTS:\n{"":=>27}{"":<48}{"":=>38}', 
               f'{"1st ITERATION VALUES:":<75}2nd ITERATION VALUES:\n{"":=<46}{"":<29}{"":=<46}']
    f17_spec = [f"{'Roof (Truss) Span:':>27}{_span:>3} m", f"{'Truss Spacing:':>27}{_S_t:>3} m",
                f"{'Thickness of Concrete Slab:':>27}{_T_s:>6} m", f"{'Snow Accumulation:':>27}{_h_s:>5} m"]
    constants = [f"{'Unit Weight of Snow (NBCC):':>38}{_ɣ_s:>5} kN/m³", 
                 f"{'Additional Weight of Wet Snow (NBCC):':>38}{_S_r:>5} kN/m²", 
                 f"{'Weight of the Steel Deck:':>38}{_W_sd:>5} kN/m²",
                 f"{'Steel Deck Flute Height:':>38}{_h_d:>7} m",
                 f"{'Number of Trusses supporting the roof:':>38}{_ROOF_SUPPORTS:>3}",
                 f"{'Unit Weight of Concrete:':>38}{_ɣ_c:>3} kN/m³",
                 f"{'Weight of Built-up Roof:':>38}{_W_r:>6} kN/m²",
                 f"{'Joist Spacing:':>38}{_S_j:>5} m"]
    first_iter = [f"{'Snow Load / Unit Area Roof:':>46}{f'{_S_s:.2f}':>7} kN/m²",
          f"{'Total Snow Load:':>46}{f'{_S:.2f}':>7} kN/m²",
          f"{'Roof Length*:':>46}{_ROOF_LENGTH:>4} m",
          f"{'Weight of Concrete Cover:':>46}{f'{_W_cc:.3f}':>8} kN/m²",
          f"{'Weight of Concrete between S.D. Flutes:':>46}{f'{_W_cf:.3f}':>8} kN/m²",
          f"{'TOTAL Weight of Roof:':>46}{f'{_W_T:.3f}':>8} kN/m²",
          f"{'Total Factored Load:':>46}{f'{_w_f:.4f}':>9} kN/m²",
          f"{'Factored, Uniform Distributed Load (Joists):':>46}{f'{_UDL_fj:.5f}':>10} kN/m",
          f"{'Factored Support Reaction (each end of joist):':>46}{f'{_P_fj:.6f}':>11} kN",
          f"{'Factored Load at each joint on the main Truss:':>46}{f'{_P_f:.5f}':>10} kN"]
    second_iter = [f"{'Open-web Steel Joist, self-weight Dead Load:':>46}{f'{_w_swj:.3f}':>8} kN/m",
          f"{'Total Dead Load per Open-web Steel Joist:':>46}{f'{_D_j:.3f}':>8} kN",
          f"{'Adjusted factored Load / joint on main Truss:':>46}{f'{_P_fa:.5f}':>10} kN"]
    
    print(f'{new_screen}{"":=<133}\n{"CURRENT VARIABLE VALUES":^133}\n{"":=>133}\n\n')
    print(headers[0])
    max_length = max([len(constants), len(f17_spec)])
    for i in range(max_length):
        try:
            print(left_fmt(f17_spec[i]), end="")
        except:
            print(left_fmt(), end="")
        try:
            print(constants[i])
        except:
            print()
    print(f'\n\n{f"{empty:=>84}":^133}\n{"CALCULATED VALUES":^133}\n{f"{empty:=>84}":^133}')
    print(headers[1])
    max_length = max([len(first_iter), len(second_iter)])
    for i in range(max_length):
        try:
            print(left_fmt(first_iter[i]), end="")
        except:
            print(left_fmt(), end="")
        try:
            print(second_iter[i])
        except:
            print()
    print('\n*not accounting for non-zero width of truss section members\n\n')
    if not first_run:
        while True:
            if from_update:
                user_choice = 'u'
            else:
                user_choice = input('[U]pdate variable values(s) or Return to [M]ain Menu [U / M]? ').lower()
            if user_choice == 'u' or user_choice == 'update':
                if not update_variable_values(True):
                    break
            elif user_choice == 'm' or re.findall('main', user_choice) or re.findall('menu', user_choice):
                break
            else:
                input("Invalid response. Please enter either <U> or <M> to make your selection.\n"
                      "Press <ENTER> to Continue...")
    
    

def calculate_HSS_radii():
    
    interm_vals = []

    for section in sections:
        interm_vals.append(re.findall('(?<=[HSS|X])\d*\.?\d*(?=X|$)', section))
        sections[section]['out_rad'] = (
            sections[section]['Area(mmÂ²)'] - 2 * float(interm_vals[-1][2]) * (
                float(interm_vals[-1][0]) + float(interm_vals[-1][1]) - π * float(interm_vals[-1][2]) / 2)) / (2 * float(interm_vals[-1][2]) * (π - 4))
        sections[section]['in_rad'] = sections[section]['out_rad'] - float(interm_vals[-1][2])



def imp_lookup_table_HSS(filename)->None:
    
    global sections
    sections = {}
    
    headers, data = [], []
    with open(filename) as file:
        file.readline()
        headers = file.readline().replace(" ", "").replace("\n", "").split(sep=",")
        for line in file:
            data.append(line.replace(" ", "").split(sep=","))
            data[-1][-1].replace("\n", "")
    
    for entry in data:
        sections[entry[0]] = {}
        for i in range(1, len(headers)):
            sections[entry[0]][headers[i]] = float(entry[i])
            
    for section in sections:
        sections[section]['T_r'] = _φ * _σ_y * sections[section]['Area(mmÂ²)'] / 1000
    
    calculate_HSS_radii()



def imp_lookup_table_OWSJ(filename)->None:
    global joists
    joists = {}
    
    headers, data, i, j = [], [], 0, 0
    
    with open(filename) as file:
        for k in range(2):
            file.readline()
        for k in range(4):
            headers.append(file.readline().strip().replace("\n", ""))
        for line in file:
            if j == 0:
                j += 1
                continue
            if j == 1:
                data.append([line.strip().replace("\n", "")])
                j += 1
                continue
            data[i] += [line.strip().replace("\n", "")]
            if j == 4:
                j = 0
                i += 1
                continue
            j += 1
    
    for entry in data:
        joists[data.index(entry)] = {}
        for element in entry:
            joists[data.index(entry)][headers[entry.index(element)]] = element



def calculate_total_mass(obj: '__main__.truss')->None:
    
    for member_type in obj.m_weights:
        for member in member_type:
            obj.total_weight += member
    obj.total_weight *= 2
    if len(obj.m_weights[2]):
        obj.total_weight -= obj.m_weights[2][-1]
    if obj.typ == 3:
        obj.total_weight -= obj.m_weights[0][3]


    
def calculate_member_mass(obj: '__main__.truss', member_length: list)->None:
    
    for i in range(len(obj.m_sections)):
        for j in range(len(obj.m_sections[i])):
            if type(obj.m_sections[i][j]) == list:
                obj.m_weights[i].append(sections[obj.m_sections[i][j][0]]['DeadLoad(kN/m)'] * member_length[i])
                continue
            obj.m_weights[i].append(sections[obj.m_sections[i][j]]['DeadLoad(kN/m)'] * member_length[i])
            
    calculate_total_mass(obj)



def assign_HSS_sections(obj: '__main__.truss')->None:
    
    _SINE, _COSINE = obj.m_dimensions['angle'][1], obj.m_dimensions['angle'][2]
    member_length = [0, obj.m_dimensions['Diag'], obj.m_dimensions['Vert'], obj.m_dimensions['Diag'], 0]
    if obj.typ == 3:
        member_length[0] = member_length[4] = obj.m_dimensions['Horz']
    elif obj.typ == 2:
        member_length[4] = obj.m_dimensions['Bot']
        member_length[0] = obj.m_dimensions['Sj']
    else:
        member_length[0] = member_length[4] = obj.m_dimensions['Sj']
    
    if obj.typ == 1:
        for i in range(len(member_length)):
            if i in [0, 3]:
                        
                for section in sections:
                    _σ_e = (π**2 * _E) / ((_K * member_length[i] / sections[section]['r(mm)'])**2)
                    λ = sqrt(_σ_y / _σ_e)
                    f = 1 / ((1 + λ**(2 * _n))**(1 / _n))
                    sections[section]['C_r'] = (sections[section]['T_r'] * f)
                    sections[section]['KLr'] = (_K * member_length[i] / sections[section]['r(mm)']) <= 200
                
                for k in range(len(obj.m_forces[i])):
                    if obj.m_forces[i][k]:
                        candidates, c2, c3 = {}, {}, []
                        for section in sections:
                            if obj.m_forces[i][k] > sections[section]['C_r']:
                                continue
                            """
                            if not sections[section]['KLr']:
                                continue
                            """                    
                            candidates[sections[section]['DeadLoad(kN/m)']] = section
                            c2[section] = sections[section]['DeadLoad(kN/m)']
                        op_sec = list(c2.values())
                        op_sec.sort()
                        j = 0
                        for item in c2.items():
                            if item[1] == op_sec[0]:
                                c3.append(item[0])
                        if len(c3) > 1:
                            ts1, ts2 = {}, {}
                            for entry in c3:
                                ts1[entry] = sections[entry]['C_r']
                                ts2[sections[entry]['C_r']] = entry
                            fc = list(ts1.values())
                            fc.sort()
                            if fc.count(fc[0]) > 1:
                                options = []
                                for item in ts1.items():
                                    if item[1] == fc[0]:
                                        options.append(item[0])
                                obj.m_sections[i].append(options)
                                continue                    
                            obj.m_sections[i].append(ts2[fc[-1]])
                            continue
                        if not candidates:
                            obj.total_weight = None
                            return                            
                        obj.m_sections[i].append(candidates[op_sec[0]])
                continue
            
            for k in range(len(obj.m_forces[i])):
                if obj.m_forces[i][k]:
                    candidates, c2, c3 = {}, {}, []
                    for section in sections:
                        if obj.m_forces[i][k] > sections[section]['T_r']:
                            continue
                        candidates[sections[section]['DeadLoad(kN/m)']] = section
                        c2[section] = sections[section]['DeadLoad(kN/m)']
                    op_sec = list(c2.values())
                    op_sec.sort()
                    j = 0
                    for item in c2.items():
                        if item[1] == op_sec[0]:
                            c3.append(item[0])
                    if len(c3) > 1:
                        ts1, ts2 = {}, {}
                        for entry in c3:
                            ts1[entry] = sections[entry]['T_r']
                            ts2[sections[entry]['T_r']] = entry
                        fc = list(ts1.values())
                        fc.sort()
                        if fc.count(fc[0]) > 1:
                            options = []
                            for item in ts1.items():
                                if item[1] == fc[0]:
                                    options.append(item[0])
                            obj.m_sections[i].append(options)
                            continue
                        obj.m_sections[i].append(ts2[fc[-1]])
                        continue
                    if not candidates:
                        obj.total_weight = None
                        return
                    obj.m_sections[i].append(candidates[op_sec[0]])
                    
        for section in sections:
            _σ_e = (π**2 * _E) / ((_K * member_length[2] / sections[section]['r(mm)'])**2)
            λ = sqrt(_σ_y / _σ_e)
            f = 1 / ((1 + λ**(2 * _n))**(1 / _n))
            sections[section]['C_r'] = (sections[section]['T_r'] * f)
            sections[section]['KLr'] = (_K * member_length[2] / sections[section]['r(mm)']) <= 200
        
        i = 2
        for k in range(1):
            if obj.m_forces[i][k]:
                candidates, c2, c3 = {}, {}, []
                for section in sections:
                    if obj.m_forces[i][k] > sections[section]['C_r']:
                        continue
                    """
                    if not sections[section]['KLr']:
                        continue
                    """                    
                    candidates[sections[section]['DeadLoad(kN/m)']] = section
                    c2[section] = sections[section]['DeadLoad(kN/m)']
                op_sec = list(c2.values())
                op_sec.sort()
                j = 0
                for item in c2.items():
                    if item[1] == op_sec[0]:
                        c3.append(item[0])
                if len(c3) > 1:
                    ts1, ts2 = {}, {}
                    for entry in c3:
                        ts1[entry] = sections[entry]['C_r']
                        ts2[sections[entry]['C_r']] = entry
                    fc = list(ts1.values())
                    fc.sort()
                    if fc.count(fc[0]) > 1:
                        options = []
                        for item in ts1.items():
                            if item[1] == fc[0]:
                                options.append(item[0])
                        obj.m_sections[i][0] = options
                        continue                    
                    obj.m_sections[i][0] = ts2[fc[-1]]
                    continue
                if not candidates:
                    obj.total_weight = None
                    return                            
                obj.m_sections[i][0] = candidates[op_sec[0]]       


    else:
        for i in range(len(member_length)):
            if i in [0, 2, 3]:
                        
                for section in sections:
                    _σ_e = (π**2 * _E) / ((_K * member_length[i] / sections[section]['r(mm)'])**2)
                    λ = sqrt(_σ_y / _σ_e)
                    f = 1 / ((1 + λ**(2 * _n))**(1 / _n))
                    sections[section]['C_r'] = (sections[section]['T_r'] * f)
                    sections[section]['KLr'] = (_K * member_length[i] / sections[section]['r(mm)']) <= 200
                
                for k in range(len(obj.m_forces[i])):
                    if obj.m_forces[i][k]:
                        candidates, c2, c3 = {}, {}, []
                        for section in sections:
                            if obj.m_forces[i][k] > sections[section]['C_r']:
                                continue
                            """
                            if not sections[section]['KLr']:
                                continue
                            """                    
                            candidates[sections[section]['DeadLoad(kN/m)']] = section
                            c2[section] = sections[section]['DeadLoad(kN/m)']
                        op_sec = list(c2.values())
                        op_sec.sort()
                        j = 0
                        for item in c2.items():
                            if item[1] == op_sec[0]:
                                c3.append(item[0])
                        if len(c3) > 1:
                            ts1, ts2 = {}, {}
                            for entry in c3:
                                ts1[entry] = sections[entry]['C_r']
                                ts2[sections[entry]['C_r']] = entry
                            fc = list(ts1.values())
                            fc.sort()
                            if fc.count(fc[0]) > 1:
                                options = []
                                for item in ts1.items():
                                    if item[1] == fc[0]:
                                        options.append(item[0])
                                obj.m_sections[i].append(options)
                                continue                    
                            obj.m_sections[i].append(ts2[fc[-1]])
                            continue
                        if not candidates:
                            obj.total_weight = None
                            return                            
                        obj.m_sections[i].append(candidates[op_sec[0]])
                continue
    
            for k in range(len(obj.m_forces[i])):
                if obj.m_forces[i][k]:
                    candidates, c2, c3 = {}, {}, []
                    for section in sections:
                        if obj.m_forces[i][k] > sections[section]['T_r']:
                            continue
                        candidates[sections[section]['DeadLoad(kN/m)']] = section
                        c2[section] = sections[section]['DeadLoad(kN/m)']
                    op_sec = list(c2.values())
                    op_sec.sort()
                    j = 0
                    for item in c2.items():
                        if item[1] == op_sec[0]:
                            c3.append(item[0])
                    if len(c3) > 1:
                        ts1, ts2 = {}, {}
                        for entry in c3:
                            ts1[entry] = sections[entry]['T_r']
                            ts2[sections[entry]['T_r']] = entry
                        fc = list(ts1.values())
                        fc.sort()
                        if fc.count(fc[0]) > 1:
                            options = []
                            for item in ts1.items():
                                if item[1] == fc[0]:
                                    options.append(item[0])
                            obj.m_sections[i].append(options)
                            continue
                        obj.m_sections[i].append(ts2[fc[-1]])
                        continue
                    if not candidates:
                        obj.total_weight = None
                        return
                    obj.m_sections[i].append(candidates[op_sec[0]])
    
    calculate_member_mass(obj, member_length)



def fast_calculate_member_forces(no_joists, angle, **kwargs)->list:
    
    single = kwargs.get("single")
    index = []
    
    if single:
        truss_type = kwargs.get("typ")        
        return_val = [truss(no_joists, angle, truss_type)]
        index = [0, 0, 0, 0]
    else:
        return_val = [truss(no_joists, angle, 0),
                      truss(no_joists, angle, 1),
                      truss(no_joists, angle, 2),
                      truss(no_joists, angle, 3)]
        index = [0, 1, 2, 3]
    
    _COSINE, _TANGENT = return_val[0].m_dimensions['angle'][2], return_val[0].m_dimensions['angle'][3]
    determinant = no_joists // 2
    
    #COMMON
    common_Horizontals = 0.5 * _P_f * _TANGENT
    common_Diagonals = _P_f / (2 * _COSINE)
    
    common_top_bot = []
    for i in range(determinant + 1):
        common_top_bot.append((no_joists - i) * (i + 1) * common_Horizontals)
    
    common_diags , common_verts = [], []
    for i in range(no_joists, 0, -2):
        common_diags.append(i * common_Diagonals)
        common_verts.append((i / 2) * _P_f)
    
    if not single or truss_type == 0:
    #HOWE TRUSS
        return_val[index[0]].m_forces[0] = common_top_bot * 1
        return_val[index[0]].m_forces[1] = common_diags * 1
        return_val[index[0]].m_forces[2] = common_verts[:-1] + [_P_f]
        return_val[index[0]].m_forces[3] = []
        return_val[index[0]].m_forces[4] = common_top_bot[:-1]
    
    if not single or truss_type == 1:
    #PRATT TRUSS
        return_val[index[1]].m_forces[0] = [common_top_bot[0]] * 2 + common_top_bot[1:-1]
        return_val[index[1]].m_forces[1] = [common_diags[0]]
        return_val[index[1]].m_forces[2] = [_P_f] + common_verts[2:] + [0]
        return_val[index[1]].m_forces[3] = common_diags[1:]
        return_val[index[1]].m_forces[4] = common_top_bot[1:]
        if no_joists == 1:
            return_val[index[1]].m_forces[0] = return_val[index[1]].m_forces[0][:-1]
            return_val[index[1]].m_forces[2] = return_val[index[1]].m_forces[2][:-1]
    
    if not single or truss_type == 2:
    #WARREN TRUSS (V)
        for i in range(0, determinant, 2):
            return_val[index[2]].m_forces[0].append(common_top_bot[i])
            return_val[index[2]].m_forces[0].append(common_top_bot[i])
            return_val[index[2]].m_forces[1].append(common_diags[i])
            return_val[index[2]].m_forces[3].append(common_diags[i + 1])
            return_val[index[2]].m_forces[2].append(_P_f)
            return_val[index[2]].m_forces[2].append(0)
            
            return_val[index[2]].m_forces[4].append(common_Horizontals * (i + 2) * (no_joists - (i + 1)))
            return_val[index[2]].m_forces[4].append(return_val[index[2]].m_forces[4][-1])
        if not determinant % 2:
            return_val[index[2]].m_forces[0].append(common_top_bot[-1])
            return_val[index[2]].m_forces[1].append(common_diags[-1])
            return_val[index[2]].m_forces[2].append(_P_f)
        else:
            return_val[index[2]].m_forces[4] = return_val[index[2]].m_forces[4][:-1]
    
    if not single or truss_type == 3:
    #WARREN TRUSS (Ø)
        Pf_R = determinant * _P_f
        common_determinant = int(determinant / 2 + 1)
        factor = 0
        for i in range(common_determinant):
            return_val[index[3]].m_forces[0].append((Pf_R * (1 + 2 * i) - 2 * _P_f * i ** 2) * _TANGENT)
        for i in range(c_round(determinant / 2)):
            factor += i
            return_val[index[3]].m_forces[4].append((2 * (i + 1) * Pf_R - 4 * factor * _P_f) * _TANGENT)
        for i in range(common_determinant):
            return_val[index[3]].m_forces[1].append((Pf_R - 2 * i * _P_f) / _COSINE)
        if ((no_joists - 1) / 2) % 2:
            return_val[index[3]].m_forces[3] = return_val[index[3]].m_forces[1]
        else:
            return_val[index[3]].m_forces[3] = return_val[index[3]].m_forces[1][:-1]
        return_val[index[3]].m_forces[2] = []
    
    if single:
        return return_val[0]
    return return_val         
    
    
 
def calculate_member_forces(no_joists, angle, **kwargs)->list:
    
    single = kwargs.get("single")
    index = []
    
    if single:
        truss_type = kwargs.get("typ")        
        return_val = [truss(no_joists, angle, truss_type)]
        index = [0, 0, 0, 0]
    else:
        return_val = [truss(no_joists, angle, 0),
                      truss(no_joists, angle, 1),
                      truss(no_joists, angle, 2),
                      truss(no_joists, angle, 3)]
        index = [0, 1, 2, 3]

    _SINE, _COSINE = return_val[0].m_dimensions['angle'][1], return_val[0].m_dimensions['angle'][2]
    determinant = [int(((no_joists - 1) / 2) % 2), no_joists // 2]
    Pf_R = no_joists * _P_f / 2
    
    if not single or truss_type == 0:
        #Howe Truss
        #SETUP
        return_val[index[0]].m_forces[2].append(Pf_R)
        return_val[index[0]].m_forces[1].append(return_val[index[0]].m_forces[2][-1] / _COSINE)
        return_val[index[0]].m_forces[0].append(return_val[index[0]].m_forces[1][-1] * _SINE)
        return_val[index[0]].m_forces[4].append(return_val[index[0]].m_forces[0][-1])
        
        for i in range(determinant[1] - 1):
            return_val[index[0]].m_forces[2].append(return_val[index[0]].m_forces[2][-1] - _P_f)
            return_val[index[0]].m_forces[1].append(return_val[index[0]].m_forces[2][-1] / _COSINE)
            return_val[index[0]].m_forces[0].append(return_val[index[0]].m_forces[1][-1] * _SINE + return_val[0].m_forces[0][-1])
            return_val[index[0]].m_forces[4].append(return_val[index[0]].m_forces[0][-1])
        
        #WRAPUP
        return_val[index[0]].m_forces[3] = []
        if no_joists == 1:
            return_val[index[0]].m_forces[2] = [_P_f]
            return_val[index[0]].m_forces[4] = []
        else:
            return_val[index[0]].m_forces[1].append((return_val[index[0]].m_forces[2][-1] - _P_f) / _COSINE)
            return_val[index[0]].m_forces[0].append(return_val[index[0]].m_forces[1][-1] * _SINE + return_val[index[0]].m_forces[0][-1])
            return_val[index[0]].m_forces[2].append(_P_f)        
        
    
    if not single or truss_type == 1:
        #Pratt Truss
        #SETUP
        return_val[index[1]].m_forces[1].append(Pf_R / _COSINE)
        return_val[index[1]].m_forces[2].append(_P_f)
        return_val[index[1]].m_forces[0].append(return_val[index[1]].m_forces[1][-1] * _SINE)
        return_val[index[1]].m_forces[0].append(return_val[index[1]].m_forces[1][-1] * _SINE)
        return_val[index[1]].m_forces[3].append(return_val[index[1]].m_forces[1][-1] - _P_f / _COSINE)
        return_val[index[1]].m_forces[4].append((return_val[index[1]].m_forces[1][-1] + return_val[index[1]].m_forces[3][-1]) * _SINE)
        
        for i in range(determinant[1] - 1):
            return_val[index[1]].m_forces[0].append(return_val[index[1]].m_forces[0][-1] + return_val[index[1]].m_forces[3][-1] * _SINE)
            return_val[index[1]].m_forces[2].append(return_val[index[1]].m_forces[3][-1] * _COSINE - _P_f)
            return_val[index[1]].m_forces[3].append(return_val[index[1]].m_forces[2][-1] / _COSINE)
            return_val[index[1]].m_forces[4].append(return_val[index[1]].m_forces[3][-1] * _SINE + return_val[index[1]].m_forces[4][-1])
            
        #WRAPUP
        if no_joists == 1:
            return_val[index[1]].m_forces[0] = return_val[index[1]].m_forces[0][:-1]
            return_val[index[1]].m_forces[3] = []
            return_val[index[1]].m_forces[4] = []
        else:
            return_val[index[1]].m_forces[2].append(0)
        
    if not single or truss_type == 2:
        #Warren(_V) Truss
        #SETUP
        return_val[index[2]].m_forces[1].append(Pf_R / _COSINE)
        return_val[index[2]].m_forces[0].append(return_val[index[2]].m_forces[1][-1] * _SINE)
        return_val[index[2]].m_forces[0].append(return_val[index[2]].m_forces[0][-1])
        vert_temp = (Pf_R - _P_f)
        return_val[index[2]].m_forces[3].append(vert_temp / _COSINE)
        return_val[index[2]].m_forces[4].append(return_val[index[2]].m_forces[0][-1] + return_val[index[2]].m_forces[3][-1] * _SINE)
        return_val[index[2]].m_forces[4].append(return_val[index[2]].m_forces[4][-1])
        return_val[index[2]].m_forces[2].append(_P_f * 1)
        return_val[index[2]].m_forces[2].append(0)
        
        for i in range(determinant[1] - 1):
            vert_temp -= _P_f
            if i % 2:
                return_val[index[2]].m_forces[3].append(vert_temp / _COSINE)
                return_val[index[2]].m_forces[4].append(return_val[index[2]].m_forces[4][-1] + _SINE * (return_val[index[2]].m_forces[1][-1] + return_val[index[2]].m_forces[3][-1]))
                return_val[index[2]].m_forces[4].append(return_val[index[2]].m_forces[4][-1])
                return_val[index[2]].m_forces[2].append(0)
            else:
                return_val[index[2]].m_forces[1].append(vert_temp / _COSINE)
                return_val[index[2]].m_forces[0].append((return_val[index[2]].m_forces[3][-1] + return_val[index[2]].m_forces[1][-1]) * _SINE + return_val[index[2]].m_forces[0][-1])
                return_val[index[2]].m_forces[0].append(return_val[index[2]].m_forces[0][-1])
                return_val[index[2]].m_forces[2].append(_P_f * 1)
        if not determinant[1] % 2:
            return_val[index[2]].m_forces[0] = return_val[index[2]].m_forces[0][:-1]
        else:
            return_val[index[2]].m_forces[4] = return_val[index[2]].m_forces[4][:-1]
        if no_joists == 1:
            return_val[index[2]].m_forces[2] = return_val[index[2]].m_forces[2][:-1]
            return_val[index[2]].m_forces[3] = []
            return_val[index[2]].m_forces[4] = []

    
    if not single or truss_type == 3:
        #Warren(_Ø) Truss
        #SETUP
        Pf_R = determinant[1] * _P_f
        return_val[index[3]].m_forces[1].append(Pf_R / _COSINE)
        return_val[index[3]].m_forces[0].append(return_val[index[3]].m_forces[1][-1] * _SINE)
        return_val[index[3]].m_forces[3].append(return_val[index[3]].m_forces[1][-1])
        return_val[index[3]].m_forces[4].append(return_val[index[3]].m_forces[0][-1] * 2)
        for i in range((int((determinant[1] - 4) / 2 + 1), int(determinant[1] / 2))[determinant[0]]):
            return_val[index[3]].m_forces[1].append(return_val[index[3]].m_forces[3][-1] - 2 * _P_f / _COSINE)
            return_val[index[3]].m_forces[0].append(return_val[index[3]].m_forces[0][-1] + (return_val[index[3]].m_forces[3][-1] + return_val[index[3]].m_forces[1][-1]) * _SINE)
            return_val[index[3]].m_forces[3].append(return_val[index[3]].m_forces[1][-1])
            return_val[index[3]].m_forces[4].append(return_val[index[3]].m_forces[4][-1] + 2 * return_val[index[3]].m_forces[1][-1] * _SINE)
        if not determinant[0]:
            return_val[index[3]].m_forces[1].append(0)
            return_val[index[3]].m_forces[0].append(return_val[index[3]].m_forces[0][-1] + (return_val[index[3]].m_forces[3][-1] + return_val[index[3]].m_forces[1][-1]) * _SINE)
        if no_joists == 1:
            return_val[index[3]].m_forces[0] = return_val[index[3]].m_forces[0][:-1]
            return_val[index[3]].m_forces[1] = return_val[index[3]].m_forces[1][:-1]
            return_val[index[3]].m_forces[3] = []
            return_val[index[3]].m_forces[4] = []
    
    if single:
        return return_val[0]
    return return_val
    
    

def initialize(re_initialize: bool=False)->None:
    
    global _span, _S_t, _T_s, _h_s, _ɣ_s, _S_r, _W_sd, _ROOF_SUPPORTS, _h_d, _ɣ_c, \
           _W_r, _S_j, _φ, _σ_y, _E, _n, _K, _angle_rng, _h_rng, _S_s, _S, \
           _ROOF_LENGTH, _W_cc, _W_cf, _W_T, _w_f, _UDL_fj, _P_fj, _P_f, _w_swj, \
           _D_j, _P_fa, number_joists    
    
    if not re_initialize:
        
        number_joists = 13
        
        _span = 35               #Roof width                 ( 35m)
        _S_t = 9                 #Truss spacing              ( 9m )
        _T_s = 0.07              #Thickness of Concrete Slab (70mm)
        _h_s = 0.8               #Greatest anticipated depth of snow (80cm)
        
        _ɣ_s = 3.2              #"(Unit?) Weight of Snow" -- NBCC (kN/m³)
        _S_r = 0.4              #Extra Weight of Snow loaded w/ rain -- NBCC (kN/m²)
        _W_sd = 0.1             #Steel Deck Weight (kN/m²)
        _ROOF_SUPPORTS = 11     #Number of Trusses supporting the roof
        _h_d = 0.038            #Steel Deck Flute Height (Thickness) (38mm)
        _ɣ_c = 24               #The "Unit Weight of Concrete" (kN/m³)
        _W_r = 0.31             #Weight of the "Built-up Roof" (kN/m²)
        _S_j = _span / 14       #Joist Spacing (FORMULA GIVEN)
        _φ = 0.9                #Some sort of material strength factor (Steel)
        _σ_y = 370              #Yield strength of Steel (MPa)
        _E = 200000             #Modulus of elasticity of Steel (MPa)
        _n = 1.34               #"A parameter for calculating compressive resistance" (thx)
                                    #For:  Hot Rolled, Fabricated Structural Sections
                                    #      Hollow Structural Sections (HSS)
                                    #Manufactured in accordance with CSA G40.20, Class C.
        _K = 1000               #"Effective Length Factor" (whatever that means)
                                #***!!!NOTE!!!***: The given K value is 1 (ONE),
                                #I HAVE CHANGED THE K-VALUE USED IN THE CALCULATIONS
                                #IN ORDER TO CONVERT L(m) INTO L(mm) TO KEEP UNITS 
                                #CONSISTENT AND CORRECT A CALCULATION PROBLEM WHERE
                                #f WOULD ALWAYS EVALUATE TO ~0.99999999 (EFFECTING
                                #NEGLIGIBLE CHANGE TO T_r WHEN CALCULATING C_r FROM IT).
        _angle_rng = [40, 50]   #The suggested range of angles between diagonal cross members and either the normal or horizontal within which we should look for an optimal solution
        _h_rng = [_span / 12, _span / 9]    #The suggested range of heights of the Truss within which we should search for an optimal solution
        
        
        #CALCULATED VALUES DEFINITIONS
        _S_s = _ɣ_s * _h_s                       #Snow load per unit area of roof (kN/m²)
        _S = _S_s + _S_r                        #TOTAL Snow load (incl. rain) (kN/m²)
        _ROOF_LENGTH = _ROOF_SUPPORTS * _S_t     #DOESN'T ACCT 4 THICKNESS OF TRUSSES (m)
        _W_cc = (_T_s - _h_d) * _ɣ_c             #Weight of "Concrete Cover" (kN/m²)
        _W_cf = (_h_d / 2) * _ɣ_c               #Wt. of C.C. Between S.D. Flutes (kN/m²)
        _W_T = _W_sd + _W_cc + _W_cf + _W_r     #TOTAL WEIGHT OF ROOF (kN/m²)
        _w_f = 1.25 * _W_T + 1.5 * _S           #FCTRD LOAD _W_T=Dead~, _S=snow~ (kN/m²)
        _UDL_fj = _w_f * _S_j                   #Unifrm Dist. Load (kN/m)
        _P_fj = _UDL_fj * _S_t / 2               #Fctr Sp Rctn = UDL * T-Spacing / 2 (kN)
        _P_f = 2 * _P_fj                        #Fctr Ld / Jnt on Trss (2 * Joists) (kN)
        
        _w_swj = 0.241                          #OWSJ Self-Weight Dead Load (kN/m)
        _D_j = _w_swj * _S_t                     #Total Dead Load per OWSJ
        _P_fa = _P_f + _D_j                     #Adjusted Factored load per joint on Truss, 
                                                #accounting for OWSJ self-weight DL (GIVEN IN ASSIGNMENT)
        """
        ***NOTE: kN/m² = kPa
                 DEAD LOAD = _w_T (not accounting for self-weight of truss / roof 
                                   dimension adjustment based on truss member section 
                                   thickness / self-weight open-web steel joists)
        """           
        
    else:
        _S_j = _span / 14
        _S_s = _ɣ_s * _h_s
        _S = _S_s + _S_r  
        _W_cc = (_T_s - _h_d) * _ɣ_c             
        _W_cf = (_h_d / 2) * _ɣ_c               
        _W_T = _W_sd + _W_cc + _W_cf + _W_r     
        _w_f = 1.25 * _W_T + 1.5 * _S           
        _UDL_fj = _w_f * _S_j                   
        _P_fj = _UDL_fj * _S_t / 2 
        _P_f = 2 * _P_fj
        _ROOF_LENGTH = _ROOF_SUPPORTS * _S_t
        
        _w_swj = 0.241
        _D_j = _w_swj * _S_t
        _P_fa = _P_f + _D_j        
        



def optimize(**kwargs)->None: 
    global weights_by_type, angle_matrix, optima, optimal_truss_design, \
           rank_order, optima_weights, t_count
    t0 = pctime()    
    
    display_graphs = kwargs.get('disp_graphs')
    if display_graphs == None:
        display_graphs = True
    basic_algorithm = kwargs.get('basic')
    if basic_algorithm == None:
        basic_algorithm = False
    display_result = kwargs.get('disp_res')
    if display_result == None:
        display_result = True
    
    weights_by_type, angle_matrix = [[], [], [], []], []
    #t_count = 0 --> moved to main()

    angles = []
    for angle in range(22, 69):
        angles.append(angle)
        angle *= pi / 180
        angle = [angle, sin(angle), cos(angle), tan(angle)]
        
        angle_matrix.append((fast_calculate_member_forces(number_joists, angle), 
                             calculate_member_forces(number_joists, angle))[basic_algorithm])
        
        t_count += 4
        
        for i in range(len(angle_matrix[-1])):
            assign_HSS_sections(angle_matrix[-1][i])
            weights_by_type[i].append(angle_matrix[-1][i].total_weight)
            
    fig, minima, range_limit = [], [], []
    for i in range(len(weights_by_type)):
        range_limit.append([])
        if not i:
            angles = [angles]
        for j in range(len(weights_by_type[i])):
            if weights_by_type[i][j] != None:
                range_limit[i].append(j)
                try:
                    range_limit[i].append(weights_by_type[i].index(None, j + 1))
                except:
                    range_limit[i].append(len(weights_by_type[i]))
                break
            else:
                range_limit[i] = [0, 0]
        weights_by_type[i] = weights_by_type[i][range_limit[i][0]:range_limit[i][1]]
        while len(angles) < i + 3:
            angles.append([])
        angles[i + 2] = angles[0][range_limit[i][0]:range_limit[i][1]]
        
        minima.append(angles[i + 2][weights_by_type[i].index(min(weights_by_type[i]))])
        
        if display_graphs:        
            fig.append(plt.figure())
            plt.xlabel('Angle to the Normal (degrees)')
            plt.ylabel('Weight of Truss (kN)')
            plt.title(('Howe Truss', 'Pratt Truss', 'Warren Truss (with vertical members)', 
                       'Warren Truss (without vertical members)')[i])
            plt.scatter(angles[i + 2], weights_by_type[i])
            
            if i == len(weights_by_type) - 1:
                fig.append(plt.figure())
                for j in range(len(weights_by_type)):
                    plt.xlabel('Angle to the Normal (degrees)')
                    plt.ylabel('Weight of Truss (kN)')
                    plt.title('Truss Weight vs Normal Angle of Diagonal Members by Truss Type')
                    plt.scatter(angles[j + 2], weights_by_type[j])
                plt.legend(['Howe Truss', 'Pratt Truss', 'Warren Truss (with vertical members)', 
                            'Warren Truss (without vertical members)'])        
    
    if display_graphs:
        rank_order, minima_w = [], []
        for i in range(len(minima)):
            rank_order.append(min(weights_by_type[i]))
            minima_w.append(rank_order[-1])
        rank_order.sort()
        rank_order = [rank_order.index(minima_w[0]) + 1, rank_order.index(minima_w[1]) + 1, 
                      rank_order.index(minima_w[2]) + 1, rank_order.index(minima_w[3]) + 1]
        figure_savenames = [f'{rank_order[0]}_Howe.png', f'{rank_order[1]}_Pratt.png', 
                            f'{rank_order[2]}_Warren_wV.png', f'{rank_order[3]}_Warren_NV.png', 
                            'ComparativeResults.png']
        for i in range(len(fig)):
            fig[i].savefig("C:\\Users\\shisc\\My Drive\\Personal\\2 - Carleton University\\"
                           "MENG - Year 1 - 2022-23\\W23\\ECOR 1046 F - Mechanics\\"
                          f"F17\\Resources\\CODE\\TestRecords\\{figure_savenames[i]}")
    
    fig, optima, plot_data = [], [], []
    for i in range(4):
        for j in range(3):
            determinant = False
            rng = []
            if j:
                if minima[i][j - 1] < minima[i][j]:
                    determinant = True
                if j == 1:
                    if determinant:
                        rng = [int(minima[i][j] * 100), int((minima[i][j] + 0.01) * 100)]
                    else:
                        rng = [(int(minima[i][j - 1] - 0.01) * 100), int(minima[i][j - 1] * 100)]
                else:
                    if determinant:
                        rng = [int(minima[i][j] * 1000), int((minima[i][j] + 0.001) * 1000)]
                    else:
                        rng = [int((minima[i][j - 1] - 0.001) * 1000), int(minima[i][j - 1] * 1000)]
            else:
                rng = [int((minima[i] - 1) * 10), int((minima[i] + 1) * 10)]
            mass_matrix = []
            x_vals = []
            results = []
            
            for angle in range(rng[0], rng[1]):
                angle /= 10 ** (j + 1)
                x_vals.append(angle)
                angle *= pi / 180
                angle = [angle, sin(angle), cos(angle), tan(angle)]
            
                results.append((fast_calculate_member_forces(number_joists, angle, single=True, typ=i), 
                                calculate_member_forces(number_joists, angle, single=True, typ=i))[basic_algorithm])
            
            for obj in results:
                assign_HSS_sections(obj)
                mass_matrix.append(obj.total_weight)
                
            t_count += len(results)
            
            if not j:
                minima[i] = [minima[i]]
            if j == 2:
                optima.append(results[mass_matrix.index(min(mass_matrix))])
                break                
            minima[i].append(x_vals[mass_matrix.index(min(mass_matrix))])
            
    elapsed_time = pctime() - t0
            
    optima_weights = []
    for obj in optima:
        optima_weights.append(obj.total_weight)
    
    optimal_truss_design = optima[optima_weights.index(min(optima_weights))]
    
    print(f'\n\nTrusses Solved / Evaluated: {t_count}\nTime Elapsed: {elapsed_time:.3f} seconds\n\n')
    if display_result:
        display_truss(optimal_truss_design)



def update_truss_matrix(*args)->None:
    
    if args:
        if len(args) > 1:
            raise ValueError("Function takes exactly 1 argument.")
        if type(args[0]) == list:
            if len(args[0][0]) != 4:
                raise ValueError("Function is expecting obj of format: <angle_matrix> "
                                 "or equivalent.")
            if not isinstance(args[0][0][0], type(angle_matrix[0][0])):
                raise ValueError("Function is expecting obj of format: <angle_matrix> "
                                 "or equivalent.")
            if type(args[0][0]) != list:
                raise ValueError("Function is expecting obj of format: <angle_matrix> "
                                 "or equivalent.")
            try:
                for angle_entry in range(len(args[0])):
                    for truss_entry in range(len(args[0][angle_entry])):
                        truss_matrix[truss_entry][angle_entry + 22] = args[0][angle_entry][truss_entry]
            except:
                raise RuntimeError("Unhandled Exception: update_truss_matrix() function "
                                   "failed to update from matrix argument.")
        elif isinstance(args[0], type(angle_matrix[0][0])):
            truss_matrix[args[0].typ][args[0].m_dimensions['angle'] * 180 / pi ] = args[0]
            
        else:
            raise ValueError("Function only accepts an object of type <class '__main__.truss'> "
                             "or an object of type list as arguments.")
    else:
        for angle_entry in range(len(angle_matrix)):
            for truss_entry in range(len(angle_matrix[angle_entry])):
                truss_matrix[truss_entry][angle_entry + 22] = angle_matrix[angle_entry][truss_entry]



def settings_menu()->None:
    
    global basic_force_algorithm, display_graphs, display_optimal, number_joists
    
    while True:
        print(f"{new_screen}{f'{empty:=<84}':^133}\n{'SETTINGS':^133}\n{f'{empty:=<84}':^133}\n\n\n\n")
        print("the following settings are available for configuration at this time:\n\n"
              f"\t1:\tForce Calculation Algorithm:\t\t\t\t\t\t{('Fast', 'Hard-Coded')[basic_force_algorithm]}\n\n"
              f"\t2:\tDisplay (total weight vs angle) graphs when solving variable \n\t\tset for optimal truss designs?\t\t\t\t\t\t{display_graphs}\n\n"
              f"\t3:\tDisplay the specifications of the optimal truss design for a given \n\t\tvariable set immediately upon completion of an optimization?\t\t{display_optimal}\n\n"
              f"\t4:\tNumber of joists in roof span:\t\t\t\t\t\t{number_joists}\n\n"
              "\t5:\tReturn to Main Menu\n")        
        user_choice = input("\n\nEnter the number corresponding with your selection from the menu above: ").lower()
        if user_choice == '1':
            basic_force_algorithm = not basic_force_algorithm
        elif user_choice == '2':
            display_graphs = not display_graphs
        elif user_choice == '3':
            display_optimal = not display_optimal
        elif user_choice == '4':
            user_choice = input("\nEnter number of joists to optimize for: ")
            if user_choice.isdigit() and int(user_choice) != number_joists:
                number_joists = int(user_choice)
                initialize(True)
            else:
                input(f"\n{user_choice} is not a valid response. Please enter an integer number.\n"
                      "Press <ENTER> to Continue...")                
        elif user_choice == '5' or user_choice == 'q' or user_choice == 'x' or \
             user_choice == 'exit' or user_choice == 'quit' or user_choice == 'm' \
             or user_choice == 'main' or user_choice == 'menu' or \
             user_choice.replace(" ", "") == 'mainmenu':
            break
        else:
            input(f"\n{user_choice} is not a valid response. Please enter a number from the menu.\n"
                  "Press <ENTER> to Continue...")



def update_variable_values(from_disp:bool=False)->bool:
    
    global _span, _S_t, _T_s, _h_s, _ɣ_s, _S_r, _W_sd, _ROOF_SUPPORTS, _h_d, _ɣ_c, \
               _W_r, _S_j, _φ, _σ_y, _E, _n, _K, _angle_rng, _h_rng, _S_s, _S, \
               _ROOF_LENGTH, _W_cc, _W_cf, _W_T, _w_f, _UDL_fj, _P_fj, _P_f, _w_swj, \
               _D_j, _P_fa, number_joists
    
    
    def left_fmt(*args)->str:
        
        if len(args) > 2:
            raise ValueError("Function left_fmt() takes a maximum of 2 arguments.")
        string, space = '', 75
        if len(args) == 2:
            if type(args[0]) == str and type(args[1]) == int:
                string, space = args[0], args[1]
            elif type(args[0]) == int and type(args[1]) == str:
                string, space = args[1], args[0]
            else:
                raise ValueError("Function left_fmt() can take only one argument of "
                                 "type string and one argument of type int.")
        if len(args) == 1:
            if type(args[0]) == str:
                string = args[0]
            elif type(args[0]) == int:
                space = args[0]
            else:
                raise ValueError("Function left_fmt() can take only arguments of "
                                 "type str or type int.")
        space -= len(string)
        space *= ' '
        return f"{string}{space}"    
    
    
    def where_to_go(loc: str)->str:
        
        while True:
            user_choice = input(f"\n\nYou may [U]pdate another variable value within {loc}, Return to the "
                                "full [L]isting of variable values\n "
                                "or, at any time, you may Return to the [M]ain Menu: ").lower()
            if user_choice == 'm' or re.findall('main', user_choice) or re.findall('menu', user_choice):
                return 'm'
            if user_choice == 'u' or re.findall('update', user_choice) or re.findall('change', user_choice):
                return 'u'
            if user_choice == 'l' or re.findall('list', user_choice) or re.findall('display', user_choice):
                return 'l'
            else:
                input(f"{user_choice} is an invalid response. Please enter from [1] to [4] (inclusive) "
                      "or [M] to return to the main menu\nPress <ENTER> to Continue...")
    
    
    def valid_number(usr_entry:str):
        
        if re.findall('\d*\.\d*', usr_entry):
            temp = usr_entry[:usr_entry.index('.')]
            temp += usr_entry[(usr_entry.index('.') + 1):]
            if temp.isdigit():
                return float(usr_entry)
            else:
                input(f"\n\n{user_choice} is an invalid response. Please enter from [1] to [4] (inclusive) "
                      "or [M] to return to the main menu\nPress <ENTER> to Continue...")
        elif usr_entry.isdigit():
            return int(usr_entry)
        else:
            input(f"\n\n{user_choice} is an invalid response. Please enter from [1] to [4] (inclusive) or "
                  "[M] to return to the main menu\nPress <ENTER> to Continue...")
            return usr_entry
    
    
    if not from_disp:
        full_list(f_update=True)
        
    else:
        headers = [f'{"F17 SPECIFIC VARIABLES:":<75}CONSTANTS:\n{"":=>27}{"":<48}{"":=>38}', 
                   f'{"1st ITERATION VALUES:":<75}2nd ITERATION VALUES:\n{"":=<46}{"":<29}{"":=<46}']
        while True:
            user_choice = input("\nPlease enter [1] if you would like to alter a "
                            "GROUP-SPECIFIC variable.\n(Listed as \"F17 SPECIFIC VARIABLES\", above)\n"
                            "Enter [2] if you would like to alter a \"CONSTANT\" value.\n"
                            "At any time you may Return to the [M]ain Menu: ").lower()
            if user_choice == 'm' or re.findall('main', user_choice) or re.findall('menu', user_choice):
                return False
            if user_choice == '1':
                while True:
                    f17_spec = [f"1:\t{'Roof (Truss) Span:':>27}{_span:>3} m", f"2:\t{'Truss Spacing:':>27}{_S_t:>3} m",
                                f"3:\t{'Thickness of Concrete Slab:':>27}{_T_s:>6} m", f"4:\t{'Snow Accumulation:':>27}{_h_s:>5} m"]                    
                    print(f'{new_screen}\n\n{"":=<133}\n{"CURRENT VARIABLE VALUES":^133}\n{"":=>133}\n\n')
                    print(f'{"F17 SPECIFIC VARIABLES:":<75}\n{"":=>27}')
                    for i in range(len(f17_spec)):
                        print(left_fmt(f17_spec[i]))
                    user_choice = input("\n\nEnter the number corresponding with the variable value you would like to change\n"
                                        " from the menu above, or, you may Return to the [M]ain Menu at any time: ").lower()
                    if user_choice == 'm' or re.findall('main', user_choice) or re.findall('menu', user_choice):
                        return False
                    if user_choice == '1':
                        user_choice = input(f"\nCurrent value of Roof Span is: {_span} m.\nEnter new value for Roof Span (meters): ")
                        user_choice = valid_number(user_choice)
                        if type(user_choice) == int or type(user_choice) == float:
                            if _span != user_choice:
                                _span = user_choice
                                initialize(True)                        
                            user_choice = where_to_go('F17 SPECIFIC VARIABLES')
                            if user_choice == 'u':
                                continue
                            elif user_choice == 'l':
                                full_list(f_update=True)
                                return False
                            elif user_choice == 'm':
                                return False
                        else:
                            continue
                    elif user_choice == '2':
                        user_choice = input(f"\nCurrent value of Truss Spacing is: {_S_t} m.\nEnter new value for Truss "
                                            "Spacing (meters): ")
                        user_choice = valid_number(user_choice)
                        if type(user_choice) == int or type(user_choice) == float:
                            if _S_t != user_choice:
                                _S_t = user_choice
                                initialize(True)                        
                            user_choice = where_to_go('F17 SPECIFIC VARIABLES')
                            if user_choice == 'u':
                                continue
                            elif user_choice == 'l':
                                full_list(f_update=True)
                                return False
                            elif user_choice == 'm':
                                return False
                        else:
                            continue                        
                    elif user_choice == '3':
                        user_choice = input(f"\nCurrent value of Concrete Slab Thickness is: {_T_s} m.\nEnter new value "
                                            "for Concrete Slab Thickness (meters): ")
                        user_choice = valid_number(user_choice)
                        if type(user_choice) == int or type(user_choice) == float:
                            if _T_s != user_choice:
                                _T_s = user_choice
                                initialize(True)                        
                            user_choice = where_to_go('F17 SPECIFIC VARIABLES')
                            if user_choice == 'u':
                                continue
                            elif user_choice == 'l':
                                full_list(f_update=True)
                                return False
                            elif user_choice == 'm':
                                return False
                        else:
                            continue                        
                    elif user_choice == '4':
                        user_choice = input(f"\nCurrent value of Snow Accumulation is: {_h_s} m.\nEnter new value for "
                                            "Snow Accumulation (meters): ")
                        user_choice = valid_number(user_choice)
                        if type(user_choice) == int or type(user_choice) == float:
                            if _h_s != user_choice:
                                _h_s = user_choice
                                initialize(True)                        
                            user_choice = where_to_go('F17 SPECIFIC VARIABLES')
                            if user_choice == 'u':
                                continue
                            elif user_choice == 'l':
                                full_list(f_update=True)
                                return False
                            elif user_choice == 'm':
                                return False
                        else:
                            continue                        
                    else:
                        input(f"\n\n{user_choice} is an invalid response. Please enter from [1] to [4] (inclusive) or "
                              "[M] to return to the main menu\nPress <ENTER> to Continue...")
            elif user_choice == '2':
                while True:
                    constants = [f"1:\t{'Unit Weight of Snow (NBCC):':>38}{_ɣ_s:>5} kN/m³", 
                                 f"2:\t{'Additional Weight of Wet Snow (NBCC):':>38}{_S_r:>5} kN/m²", 
                                 f"3:\t{'Weight of the Steel Deck:':>38}{_W_sd:>5} kN/m²",
                                 f"4:\t{'Steel Deck Flute Height:':>38}{_h_d:>7} m",
                                 f"5:\t{'Number of Trusses supporting the roof:':>38}{_ROOF_SUPPORTS:>3}",
                                 f"6:\t{'Unit Weight of Concrete:':>38}{_ɣ_c:>3} kN/m³",
                                 f"7:\t{'Weight of Built-up Roof:':>38}{_W_r:>6} kN/m²",
                                 f"8:\t{'Joist Spacing:':>38}{_S_j:>5} m"]                    
                    print(f'{new_screen}\n\n{"":=<133}\n{"CURRENT VARIABLE VALUES":^133}\n{"":=>133}\n\n')            
                    print(f'{"CONSTANTS:":<75}\n{"":=>38}')
                    for i in range(len(constants)):
                        print(left_fmt(constants[i]))
                    user_choice = input("\n\nEnter the number corresponding with the variable value you would like to change\n"
                                        " from the menu above, or, you may Return to the [M]ain Menu at any time: ").lower()
                    if user_choice == 'm' or re.findall('main', user_choice) or re.findall('menu', user_choice):
                        return False
                    if user_choice == '1':
                        user_choice = input(f"\nCurrent value of the Unit Weight of Snow (NBCC) is: {_ɣ_s} kN/m³.\n"
                                            "Enter new value for Unit Weight of Snow (kN/m³): ")
                        user_choice = valid_number(user_choice)
                        if type(user_choice) == int or type(user_choice) == float:                        
                            if _ɣ_s != user_choice:
                                _ɣ_s = user_choice
                                initialize(True)
                            user_choice = where_to_go('CONSTANTS')
                            if user_choice == 'u':
                                continue
                            elif user_choice == 'l':
                                full_list(f_update=True)
                                return False
                            elif user_choice == 'm':
                                return False
                        else:
                            continue                            
                    elif user_choice == '2':
                        user_choice = input(f"\nCurrent value of Additional Weight of WET Snow is: {_S_r} kN/m².\nEnter "
                                            "new value for Additional Weight of WET Snow (kN/m²): ")
                        user_choice = valid_number(user_choice)
                        if type(user_choice) == int or type(user_choice) == float:
                            if _S_r != user_choice:
                                _S_r = user_choice
                                initialize(True)
                            user_choice = where_to_go('CONSTANTS')
                            if user_choice == 'u':
                                continue
                            elif user_choice == 'l':
                                return False
                            elif user_choice == 'm':
                                return False
                        else:
                            continue                            
                    elif user_choice == '3':
                        useer_choice = input(f"\nCurrent value of the Weight of the Steel Deck is: {_W_sd} m.\nEnter "
                                             "new value for the Weight of the Steel Deck (meters): ")
                        user_choice = valid_number(user_choice)
                        if type(user_choice) == int or type(user_choice) == float:                        
                            if _W_sd != user_choice:
                                _W_sd = user_choice
                                initialize(True)
                            user_choice = where_to_go('CONSTANTS')
                            if user_choice == 'u':
                                continue
                            elif user_choice == 'l':
                                full_list(f_update=True)
                                return False
                            elif user_choice == 'm':
                                return False
                        else:
                            continue                        
                    elif user_choice == '4':
                        user_choice = input(f"\nCurrent value of the Steel Deck Flute Heightc(/depth) is: {_h_d} m.\n"
                                            "Enter new value for the Steel Deck Flute Heightc(/depth) (meters): ")
                        user_choice = valid_number(user_choice)
                        if type(user_choice) == int or type(user_choice) == float:                        
                            if _h_d != user_choice:
                                _h_d = user_choice
                                initialize(True)
                            user_choice = where_to_go('CONSTANTS')
                            if user_choice == 'u':
                                continue
                            elif user_choice == 'l':
                                full_list(f_update=True)
                                return False
                            elif user_choice == 'm':
                                return False
                        else:
                            continue                            
                    elif user_choice == '5':
                        input(f"\nCurrent value of the Number of Trusses Supporting the Roof is: {_ROOF_SUPPORTS}.\n"
                              "However, while this value is noted and displayed, it has NO EFFECT on ANY truss calculations."
                              " Thus it is not a variable whose value can be adjusted.\n"
                              "Press <ENTER> to Continue...")
                    elif user_choice == '6':
                        user_choice = input(f"\nCurrent value of the Unit Weight of Concrete is: {_ɣ_c} kN/m³.\nEnter new "
                                            "value for the Unit Weight of Concrete (kN/m³): ")
                        user_choice = valid_number(user_choice)
                        if type(user_choice) == int or type(user_choice) == float:                        
                            if _ɣ_c != user_choice:
                                _ɣ_c = user_choice
                                initialize(True)
                            user_choice = where_to_go('CONSTANTS')
                            if user_choice == 'u':
                                continue
                            elif user_choice == 'l':
                                full_list(f_update=True)
                                return False
                            elif user_choice == 'm':
                                return False
                        else:
                            continue                            
                    elif user_choice == '7':
                        user_choice = input(f"\nCurrent value of Roof Span is: {_W_r}.\nEnter new value for Roof Span (meters): ")
                        user_choice = valid_number(user_choice)
                        if type(user_choice) == int or type(user_choice) == float:                        
                            if _W_r != user_choice:
                                W_r = user_choice
                                initialize(True)
                            user_choice = where_to_go('CONSTANTS')
                            if user_choice == 'u':
                                continue
                            elif user_choice == 'l':
                                full_list(f_update=True)
                                return False
                            elif user_choice == 'm':
                                return False
                        else:
                            continue
                    elif user_choice == '8':
                        input(f"\nCurrent value of the Joist Spacing is: {_S_j}.\n"
                              "However, because this value is technically a calculated value (based on the span and the number of "
                              "joists) and adjusting it independently of those other variables\n"
                              "would require the deep consideration of its effects on every area of truss calculation, this feature"
                              " has been omitted for now.\n"
                              "It may, however, be available at some point in the future...\n"
                              "Press <ENTER> to Continue...")
                    else:
                        input(f"\n\n{user_choice} is an invalid response. Please enter from [1] to [8] (inclusive) or [M]\n"
                              "Press <ENTER> to Continue...")                
            else:
                input(f"\n\n{user_choice} is an invalid response. Please enter [1], [2] or [M]\n"
                      "Press <ENTER> to Continue...")
                full_list(f_update=True)
                return False



def main_menu()->bool:
    
    global disp_from_disp, t_count
    disp_from_disp = False
    print(f"{new_screen}{f'{empty:=<84}':^133}\n{'MAIN MENU':^133}\n{f'{empty:=<84}':^133}")
    print("\n\n\t1:\tOptimize\n\t2:\tRun Diagnostics\n\t3:\tDisplay Current Variable Values\n"
          "\t4:\tDisplay Optimal Truss\n\t5:\tDisplay Truss \\\n\t\tCalculate single truss\n\t6:\tUpdate "
          "Variable Values\n\t7:\tRestore Default Variable Values\n\t8:\t"
          "Settings\n\t9:\t(Q)uit \ (E)xit\n")
    user_choice = input("Please enter the number corresponding with your choice from the menu above: ").lower()
    if user_choice == '1':
        user_choice = input(("\n\n\nRun Optimization [Y / N]? ", 
                             "\n\n\nOptimization has already been run. Overwrite existing truss data [Y / N]? ")[bool(truss_matrix[0])]).lower()
        t_count = 0
        if truss_matrix[0]:
            if user_choice == 'n' or user_choice == 'no':
                user_choice = input("\nRun Optimization without storing resultant data [Y / N]? ").lower()
                if user_choice == 'y' or 'yes':
                    optimize(disp_res=display_optimal, disp_graphs=display_graphs)
                    input("Press <ENTER> to return to main menu.")
            else:
                optimize(disp_res=display_optimal, disp_graphs=display_graphs)
                update_truss_matrix()
                input("Press <ENTER> to return to main menu.")
        else:
            if user_choice == 'y' or user_choice == 'yes':
                optimize(disp_res=display_optimal, disp_graphs=display_graphs)
                update_truss_matrix()
                input("Press <ENTER> to return to main menu.")
    elif user_choice == '2':
        user_choice = input(("\n\n\nRun Diagnostics on Force Calculation Algorithm [Y / N]? ", 
                             "\n\n\nOptimization has already been run. Overwrite existing truss data [Y / N]? ")[bool(truss_matrix[0])]).lower()
        if truss_matrix[0]:
            if user_choice == 'y' or user_choice == 'yes':
                user_choice = input("\nRun Diagnostics without storing resultant optimization data [Y / N]? ").lower()
                if user_choice == 'y' or 'yes':
                    diagnostics.forces(no_update=True)
            else:
                diagnostics.forces()
        if user_choice == 'y' or user_choice == 'yes':
            diagnostics.forces()
    elif user_choice == '3':
        full_list()
    elif user_choice == '4':
        try:
            disp_from_disp = True
            display_truss(optimal_truss_design)
        except:
            input("\n\nNo Optimization results loaded into memory. Run optimization or load saved results from file.\n"
                  "Press <ENTER> to Continue...")
    elif user_choice == '5':
        if truss_matrix[0]:
            disp_from_disp = True
            disp_truss_menu()
        else:
            input("\n\nNo Optimization results loaded into memory. Run optimization or load saved results from file.\n"
                  "Press <ENTER> to Continue...")            
    elif user_choice == '6':
        update_variable_values()
    elif user_choice == '7':
        user_choice = input("\n\nAre you SURE? Restore ALL Variable Values to Defaults [Y / N]?  ").lower()
        if user_choice == 'y' or user_choice == 'yes':
            initialize()
            print(f"{new_screen}")
            for i in range(3):
                print(". ", end="")
                sleep(1)
            input("Variables reinitialized to default (group F17 / ECOR1046F: Project) values.\n"
                  "Press <ENTER> to Continue...")
    elif user_choice == '8':
        settings_menu()
    elif user_choice == '9' or user_choice == 'e' or user_choice == 'q' \
         or user_choice == 'exit' or user_choice == 'quit':
        clean_exit = input("\nRetain Optimization data / program variables? (Accessible from wing console) [Y / N]?  ").lower()
        if clean_exit == 'n' or clean_exit == 'no':
            clean_exit = True
        else:
            clean_exit = False
        return True
    else:
        input(f"{user_choice} is not a valid response. Please select the number corresponding with the your "
              "choice from the menu.\n Press <ENTER> to Continue...")    





def main()->None:
    
    global t_count, t0, clean_exit, empty, diagnostics, new_screen, translate, \
           truss_matrix, basic_force_algorithm, display_graphs, display_optimal
    t_count, empty, clean_exit, new_screen, t0 = 0, '', False, "\n", 0
    truss_matrix = {0: {}, 1: {}, 2: {}, 3: {}}
    new_screen *= 43
    basic_force_algorithm, display_graphs,  display_optimal = False, True, True
    diagnostics = diag_obj()
    
    initialize()
    imp_lookup_table_HSS('HSS_lookup_table.csv')
    imp_lookup_table_OWSJ('OWSJ_lookup_table.csv')
    
    print(f"\n\n{'':*>24}\nWELCOME TO TRUSS-ty-CALC v1.2.1\n{'':*>24}\n\n\n")
    
    while True:
        full_list(first_rn=True)
        user_choice = input("\n\nAccept default variable values [Y / N]? ").lower()
        
        if user_choice == 'y' or user_choice == 'yes':
            break
        elif user_choice == 'n' or user_choice == 'no':
            update_variable_values()
            break
        else:
            input(f"\n{user_choice} is not a valid response. Please select <Y> or <N>.\nPress <ENTER> to Continue...")
    
    quit = False
    while not quit:
        quit = main_menu()
        
    
    
    
if __name__ == '__main__':
    
    main()
    print("\n\nThank you for using TRUSS-ty-CALC! Progam will now exit.")
    sleep(3)
    
    variables = [ t_count, t0, _span, _S_t, _T_s, _h_s, \
            _ɣ_s, _S_r, _W_sd, _ROOF_SUPPORTS, _h_d, _ɣ_c, _W_r, _S_j, _φ, _σ_y, \
            _E, _n, _K, _angle_rng, _h_rng, _S_s, _S, _ROOF_LENGTH, _W_cc, _W_cf, \
            _W_T, _w_f, _UDL_fj, _P_fj, _P_f, _w_swj, _D_j, _P_fa, joists, sections, \
            assign_HSS_sections, calculate_HSS_radii, calculate_member_forces, \
            calculate_member_mass, calculate_total_mass, c_round, display_truss, \
            fast_calculate_member_forces, full_list, imp_lookup_table_HSS, \
            imp_lookup_table_OWSJ, initialize, main, optimize, truss]
    
    if clean_exit:
        try:
            del weights_by_type, angle_matrix, optima, optimal_truss_design, \
                rank_order, optima_weights
        except:
            pass
        try:
            del convert_mag, measure_runtime
        except:
            pass
        
        for variable in variables:
            try:
                del variable
            except:
                continue
        collect()
