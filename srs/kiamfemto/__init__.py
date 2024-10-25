"""
Проект kiam-femto создан для...
"""
import sys
print(sys.path[0])
sys.path.insert(1, f"{sys.path[0]}/kiamfemto")
sys.path.insert(1, f"{sys.path[0]}/test")

from config import *
from cosmetic import *
from dynamics import *
from simulation import *
from gnc_systems import *
from interface import *
from my_math import *
from my_plot import *
from primary_info import *
from spacecrafts import *

my_print(f"Инициализация проекта kiam-femto | Контекст: {__name__}", color="g")
