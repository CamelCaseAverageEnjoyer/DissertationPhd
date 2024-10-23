import sys
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

# from cosmetic import *
from my_plot import *

class Example(QWidget):
    def __init__(self, o):
        super().__init__()
        self.o = o
        self.config_choose_n = 0
        self.path = "../images/"
        self.n = 100
        self.wb = 50
        self.name_type_func = [[['', '', None] for _ in range(self.n)] for _ in range(self.n)]
        self.initUI()

    def initUI(self):
        grid = QGridLayout()
        self.setLayout(grid)
        grid.setSpacing(10)

        # Редактирование сетки
        self.buttons_and_labels()

        positions = [(i, j) for i in range(self.n) for j in range(self.n)]
        names = [self.name_type_func[i][j][0] for i in range(self.n) for j in range(self.n)]
        types = [self.name_type_func[i][j][1] for i in range(self.n) for j in range(self.n)]
        funcs = [self.name_type_func[i][j][2] for i in range(self.n) for j in range(self.n)]
        for position, name, t, f in zip(positions, names, types, funcs):
            if name == '':
                continue
            if t == 'button':
                if 'jpg' in name or 'png' in name:
                    button = QPushButton('')
                    button.setIcon(QIcon(name))
                    button.setIconSize(QSize(self.wb, self.wb))
                else:
                    button = QPushButton(name)
                if f is not None:
                    # button.clicked.connect(f)
                    button.pressed.connect(f)
                grid.addWidget(button, *position)
            elif t == 'label':
                label = QLabel(name)
                label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
                grid.addWidget(label, *position)


        # Редактирование окна
        self.move(0, 0)
        self.resize(920, 580)
        self.setWindowTitle('kiam-femto')
        self.setStyleSheet('background-color: grey;')
        self.setWindowIcon(QIcon(self.path + "wizard.png"))
        self.show()

    def buttons_and_labels(self):
        # Ручное редактирование сетки
        self.name_type_func[0][0] = [self.path + "robot1.png", "button", talk]
        self.name_type_func[1][0] = [self.path + "integral.png", "button", self.main_run]

        # Заполнение выбора config
        self.name_type_func[0][4] = ["Сохранить текущие параметры", "button", self.o.v.save_params]
        print(self.o.v.config_choose)
        for i in range(len(self.o.v.config_choose)):
            self.name_type_func[i+1][3] = ["Параметры: " + self.o.v.config_choose.iloc[i, 0], "label", None]
            self.name_type_func[i+1][4] = ["Загрузить", "button", self.o.v.load_params]

    def main_run(self):
        self.o.integrate(t=self.o.v.TIME, animate=False)

        # Вывод результатов
        tmp = np.array([np.linalg.norm(self.o.f.line_difference[0][i])
                        for i in range(len(self.o.f.line_difference[0]))])
        print(f"Математическое ожидание: {tmp.mean()}, Среднее отклонение: {tmp.std()}")
        talk_decision(cnd=self.o.v.IF_TALK)
        # plot_all(self.o)
        # plot_sigmas(self.o)
        plot_distance(self.o)  # Можно не комментировать

def interface_window(o):
    app = QApplication(sys.argv)
    window = Example(o=o)
    return app, window
