#!/usr/bin/env python3


from math import log as ln
from pprint import pprint
from random import gauss, randrange
from getpass import getpass


def distance_modulus(distance_in_pc):
    return 5*log(distance_in_pc) - 5


def gauss_mag(mean, sigma):
    return gauss(mean, sigma)


def mag_v(N1):
    M_VSol = 24/5
    return M_VSol + 1/5*(N1 - 7/2)


def mean_distance(dist_range):
    return dist_range.start + (dist_range.stop - dist_range.start)/2


def metallicity(N2):
    return (2*N2 + 1)/12


def metallicity_dimming(N2):
    return 87/100*log(1/metallicity(N2))


def log(number, base=10): 
    return ln(number, base)


def print_table(sky, per_line=6):
    BOLD, RESET  = '\033[1m', '\033[0m'
    bold_if = lambda x, t: f'{BOLD}{x:<5.4g}{RESET}' if x < t else f'{x:<5.4g}'
    header = f'{"Region":<6} {"n":>3} {"d(pc)":>7} '
    header+= f'{"μ":>7} {"10−μ":>7}  M_V (bold = visible)'
    print(f'{header}\n{"-"*len(header)}')
    for region in sky:
        limit = 10 - sky[region]['mu']
        vals = [bold_if(m, limit) for m in sky[region]['M_V']]
        chunks = [vals[i:i + per_line] for i in range(0, len(vals), per_line)]
        line1 = f'{region:<6} {sky[region]["n"]:>3} {int(sky[region]["d"]):>7}'
        line1+= f' {sky[region]["mu"]:>7.4g} {limit:>7.4g}  '
        line2 = f'{"":<6} {"":>3} {"":>7} {"":>7} {"":>7}  '
        line2+= f'visible: {len(sky[region]["sample"])}/{sky[region]["n"]}'
        if chunks:
            print(line1 + ' '.join(chunks[0]))
            for chunk in chunks[1:]: print(' '*len(line1) + ' '.join(chunk))
        else: print(line1)
        print(line2)
    return


def uniform_density(region1, region2):
    V = [r.stop**3 - r.start**3 for r in [region1['range'], region2['range']]]
    return region1['n']*(V[1]/V[0])


if __name__ == '__main__':

    sky = {
        'A': {'range': range( 70,  90)},
        'B': {'range': range( 90, 110)},
        'C': {'range': range(110, 130)}}
    
    sky['B']['n'] = 50
    M_V_solar, sigma = 24/5, 3/10
    for region in sky:
        sky[region]['n'] = round(uniform_density(sky['B'], sky[region]))
        sky[region]['d'] = mean_distance(sky[region]['range'])
        sky[region]['mu'] = distance_modulus(sky[region]['d'])
        sky[region]['M_V'] = [
            gauss_mag(M_V_solar, sigma) for x in range(sky[region]['n'])]
        sky[region]['sample'] = [
            x for x in sky[region]['M_V'] if x < 10 - sky[region]['mu']]
    M_V_all = [mag for region in sky for mag in sky[region]['M_V']]
    M_V_sample = [mag for region in sky for mag in sky[region]['sample']]
    avg_M_V = sum(M_V_all)/len(M_V_all)
    avg_M_V_sample = sum(M_V_sample)/len(M_V_sample)
    delta_M = avg_M_V - avg_M_V_sample
    
    print('Part (c)')
    print_table(sky)
    print(f'\tΔM̄ = M̄_all − M̄_sample = {delta_M:.4g} mag')
    getpass('\nPress Enter to continue...\n')

    avg_d = sum([sky[x]['d']*sky[x]['n'] for x in sky])/len(M_V_all)
    avg_d_sample = sum(
            [sky[x]['d']*len(sky[x]['sample']) for x in sky])/len(M_V_sample)
    print(
        f'Part (d)\nd̄_all = {avg_d:.4g} pc\n'
        f'd̄_sample = {avg_d_sample:.4g} pc')
    getpass('\nPress Enter to continue...\n')

    hat_d_sample = sum([
        sky[x]['d']*10**((y - avg_M_V)/5) 
        for x in sky for y in sky[x]['sample']])/len(M_V_sample)
    delta_d = avg_d_sample - hat_d_sample

    print(
        f'Part (e)\n'
        f'If we assume our sample stars have the average luminosity for all '
        f'stars and\nthen calculated their distance from their apparent '
        f'magnitudes then \n\that d̄_sample = {hat_d_sample:.4g}.\nWe would '
        f'be {"over" if delta_d < 0 else "under"}estimating the '
        f'average distance by {abs(delta_d):.4g} parsecs')
    getpass('\nPress Enter to continue...\n')

    for region in sky:
        sky[region]['N2'] = [randrange(1, 7) 
            for x in range(len(sky[region]['M_V']))]
        sky[region]['M_V'] = [x + metallicity_dimming(y) 
            for x, y in zip(sky[region]['M_V'], sky[region]['N2']) ]
        sky[region]['sample'] = [x 
            for x in sky[region]['M_V'] if x < 10 - sky[region]['mu']]
    M_V_all = [mag for region in sky for mag in sky[region]['M_V']]
    M_V_sample = [mag for region in sky for mag in sky[region]['sample']]
    avg_M_V = sum(M_V_all)/len(M_V_all)
    avg_M_V_sample = sum(M_V_sample)/len(M_V_sample)
    delta_M = avg_M_V - avg_M_V_sample
    
    print('Part (f) - Metallicity Dimming')
    print_table(sky)
    print(f'\tΔM̄ = M̄_all − M̄_sample = {delta_M:.4g} mag')
    getpass('\nPress Enter to continue...\n')

    Z_all = [metallicity(y) for x in sky for y in sky[x]['N2']]
    Z_all = sum(Z_all)/len(Z_all)
    Z_BC_sample =[
        metallicity(sky[x]['N2'][y]) 
        for x in ['B', 'C'] for y in range(len(sky[x]['N2']))
        if sky[x]['M_V'][y] in sky[x]['sample']]
    Z_BC_sample = sum(Z_BC_sample)/len(Z_BC_sample)
    more_less = 'more' if Z_all < Z_BC_sample else 'less'

    print(f'Part (g)')
    print(f'\tZ_all: {Z_all:.4g}')
    print(f'\tZ_sample: {Z_BC_sample:.4g}')
    print(
        f'The observed outer-region sample is {more_less}\nmetal-rich '
        f'than the overall population.')


