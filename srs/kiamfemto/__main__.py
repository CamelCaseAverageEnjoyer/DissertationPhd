"""Численное моделирование космической миссии с использованием чипсатов"""
from interface import *

if __name__ == '__main__':
    o = Objects(v=Variables())

    # Интерфейс
    app, window = interface_window(o=o)
    sys.exit(app.exec_())

