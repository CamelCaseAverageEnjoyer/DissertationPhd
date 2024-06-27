"""Проверка на работоспособность установленных библиотек. Среди них:
- kiamastro
    Нет FKIAMToolbox, не может ничего запустить
- nasapy
    Есть класс Nasa, что-то делает
- celmech
    Большая библиотека. Обязательно изучить, когда займёшься уравнениями
- poliastro
    Атмосфера: https://docs.poliastro.space/en/stable/examples/Atmospheric%20models.html
- astropy
"""

from kiam_astro import trajectory, kiam
import numpy as np

# help(kiam.units('earth'))
print(kiam.units('earth')['TimeUnit'] * 24 * 2*np.pi)  # 'DistUnit', 'VelUnit', 'TimeUnit'
