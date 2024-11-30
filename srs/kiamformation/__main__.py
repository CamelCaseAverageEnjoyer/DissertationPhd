"""Численное моделирование космической миссии с использованием чипсатов"""
from interface import *
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)


if __name__ == '__main__':
    # Инициализация объектов
    o = init()
    np.set_printoptions(linewidth=300)

    # Интерфейс
    app, window = interface_window(o=o)
    sys.exit(app.exec_())

