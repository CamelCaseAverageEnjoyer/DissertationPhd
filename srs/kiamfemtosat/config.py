# >>>>>>>>>>>> Вручную настраиваемые параметры <<<<<<<<<<<<
SHAMANISM = {'KalmanQuaternionNormalize': True,   # Нормировка кватернионов в фильтре Калмана
             'KalmanSpinLimit': [True, 1e-3],  # Ограничение скорости вращения в прогнозе фильтра Калмана
             'ClohessyWiltshireC1=0': True}  # Траектории без дрейфа (зануление C1, при учёте аэродинамики поломок нет)
R_V_CubeSat_SPREAD = [0, 0]
DISTORTION = 0.
NO_LINE_FLAG = -10

# >>>>>>>>>>>> Константы <<<<<<<<<<<<
GAIN_MODES = ['isotropic', 'ellipsoid', '1 antenna', '2 antennas', '1+1 antennas']
NAVIGATIONS = ['perfect', 'near', 'random']
OPERATING_MODES = ['free_flying', 'swarm_stabilize', 'lost']
OPERATING_MODES_CHANGE = ['const', 'while_sun_visible']
MY_COLORS = ['violet', 'blueviolet', 'forestgreen', 'cornflowerblue', 'peru', 'teal', 'blueviolet', 'deeppink',
             'darksalmon', 'magenta', 'maroon', 'orchid', 'purple', 'wheat', 'tan', 'steelblue', 'forestgreen',
             'aqua', 'blue', 'beige', 'bisque', 'indigo', 'navy', 'deepskyblue', 'maroon', 'gold', 'aquamarine',
             'indigo', 'olivedrab', 'slategray', 'pink', 'salmon', 'steelblue', 'peru']
