"""
Проект kiam-femto создан для...
"""
import sys
sys.path.insert(1, f"{sys.path[0]}/kiamformation")
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
from symbolic import *

my_print(f"Инициализация проекта kiam-formation | Контекст: {__name__}", color="g")
