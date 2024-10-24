import sys
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

# from cosmetic import *
from my_plot import *

class Window(QWidget):
    def __init__(self, o):
        super().__init__()
        self.o = o
        self.path = "../images/"
        self.n = 100
        self.wb = 50
        self.textboxes = {}
        self.checkboxes = {}
        self.comboboxes = {}
        self.choices = {}
        self.name_type_func = []
        self.grid = QGridLayout()
        self.setLayout(self.grid)
        self.grid.setSpacing(10)
        self.config_choose_n = 0

        self.initUI()

    def initUI(self):
        # Очистка layout
        for i in reversed(range(self.grid.count())):
            self.grid.itemAt(i).widget().deleteLater()
        self.name_type_func = [[['', '', None, (1, 1)] for _ in range(self.n)] for _ in range(self.n)]

        # Подгрузка параметров
        self.o.v.load_params(i=self.config_choose_n)

        # Редактирование сетки
        self.buttons_and_labels()

        positions = [(i, j) for i in range(self.n) for j in range(self.n)]
        names = [self.name_type_func[i][j][0] for i in range(self.n) for j in range(self.n)]
        types = [self.name_type_func[i][j][1] for i in range(self.n) for j in range(self.n)]
        funcs = [self.name_type_func[i][j][2] for i in range(self.n) for j in range(self.n)]
        whs = [self.name_type_func[i][j][3] for i in range(self.n) for j in range(self.n)]
        for position, name, t, f, wh in zip(positions, names, types, funcs, whs):
            if name == '':
                continue
            if 'button' in t:
                if 'jpg' in name or 'png' in name:
                    button = QPushButton('')
                    button.setIcon(QIcon(name))
                    button.setIconSize(QSize(self.wb, self.wb))
                else:
                    button = QPushButton(name)
                if f is not None:
                    # button.clicked.connect(f)
                    button.pressed.connect(f)
                self.grid.addWidget(button, *position, *wh)
            elif 'label' in t:
                label = QLabel(name)
                label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
                label.setStyleSheet("border: 1px solid #C2C7CB; border-radius: 5%; padding: 10px;")
                self.grid.addWidget(label, *position, *wh)
            elif 'image' in t:
                label = QLabel()
                pix = QPixmap(name)
                label.setPixmap(pix)
                # label.setPixmap(pix.scaled(500, 500, Qt.AspectRatioMode.KeepAspectRatio,
                #                            Qt.TransformationMode.FastTransformation))
                self.grid.addWidget(label, *position, *wh)
            elif 'edit' in t:
                textbox = QLineEdit()
                textbox.setText(t.split(";")[1])
                self.textboxes[name] = textbox
                self.grid.addWidget(textbox, *position, *wh)
            elif 'check' in t:
                b = QCheckBox()
                b.setChecked(int(t.split(";")[1]) > 0)
                self.checkboxes[name] = b
                self.grid.addWidget(b, *position, *wh)
            elif 'combo' in t:
                c = QComboBox()
                c.addItems(f.split(";"))
                c.setCurrentText(t.split(";")[1])
                self.comboboxes[name] = c
                self.grid.addWidget(c, *position, *wh)

        # Изображение
        # painter = QPainter(self)
        # painter.setPen(Qt.blue)

        # n_line = 100
        # radius = 100

        rounds = self.o.v.TIME / self.o.v.MY_SEC_IN_TURN
        # trajectory = [QPoint(i*10000, i*10000) for i in range(n_line)]
        # painter.drawPolyline(QPolygon(trajectory))
        # radius * np.cos(rounds * i / n_line) + i*10,
        #                              radius * np.sin(rounds * i / n_line)

        # Редактирование окна
        self.move(0, 0)
        # self.resize(1920, 1080)
        self.setWindowTitle('kiam-femto')
        self.setStyleSheet('background-color: grey;')
        self.setWindowIcon(QIcon(self.path + "wizard.png"))
        self.show()

    def buttons_and_labels(self):
        data = self.o.v.config_choose
        params = data.iloc[self.config_choose_n, :]
        N_right = 8
        N_down = 9

        # Ручное редактирование сетки
        self.name_type_func[0][0] = [self.path + "robot1.png", "button", talk, (1, 1)]  # Поболтать
        self.name_type_func[1][0] = [self.path + "integral.png", "button", self.main_run, (1, 1)]  # Моделирование
        self.name_type_func[2][0] = [self.path + "antenna.png", "button", plot_model_gain, (1, 1)]  # Диаграммы
        self.name_type_func[3][0] = [self.path + "air.png", "button", plot_atmosphere_models, (1, 1)]
        self.name_type_func[4][0] = [self.path + "animation.png", "button", animate_reference_frames, (1, 1)]

        # Автоматическое заполнение сохранений
        self.name_type_func[0][N_right+0] = ["Параметры", "edit;Название", None, (1, 2)]
        self.name_type_func[0][N_right+2] = ["Сохранить", "button", self.local_save_params, (1, 1)]
        for i in range(len(self.o.v.config_choose)):
            self.name_type_func[i+1][N_right+0] = [data.iloc[i, 0], "label", None, (1, 1)]
            self.name_type_func[i+1][N_right+1] = [f"Загрузить", "button", lambda j=i: self.local_load_params(i=j),
                                                   (1, 1)]
            self.name_type_func[i+1][N_right+2] = [f"Удалить", "button", lambda j=i: self.local_remove_params(i=j),
                                                   (1, 1)]
        # Отрисовка
        self.name_type_func[N_down][1] = [f"Разложения сигналов: {params['MULTI_ANTENNA_SEND']} -> "
                                          f"{params['MULTI_ANTENNA_TAKE']}", "label", None, (1, 3)]
        self.name_type_func[N_down+1][1] = [self.path + "send0_take0.png", "image", None, (1, 3)]
        self.name_type_func[N_down][4] = [self.o.time_message(params['TIME']), "label", None, (1, 3)]

        # Параметры численного моделирования
        self.name_type_func[0][1] = ["Параметры численного моделирования", "label", None, (1, N_right - 2)]

        # Первый столбец
        self.name_type_func[1][1] = ["Шаг по времени dT", "label", None, (1, 1)]
        self.name_type_func[1][2] = ["dT", f"combo;{params['dT']}", "0.1;1.0;10.0", (1, 1)]
        self.name_type_func[2][1] = ["Время интегрирования T", "label", "", (1, 1)]
        self.name_type_func[2][2] = ["TIME", f"combo;{params['TIME']}", "1000.0;10000.0;100000.0", (1, 1)]
        self.name_type_func[3][1] = ["Априорная информация", "label", "", (1, 1)]
        self.name_type_func[3][2] = ["START_NAVIGATION_N",
                                     f"combo;{self.o.v.NAVIGATIONS[params['START_NAVIGATION_N']]}",
                                     ";".join(self.o.v.NAVIGATIONS), (1, 1)]
        self.name_type_func[4][1] = ["Антенны кубсатов", "label", "", (1, 1)]
        self.name_type_func[4][2] = ["GAIN_MODEL_C_N",
                                     f"combo;{self.o.v.GAIN_MODES[params['GAIN_MODEL_C_N']]}",
                                     ";".join(self.o.v.GAIN_MODES), (1, 1)]
        self.name_type_func[5][1] = ["Антенны чипсатов", "label", "", (1, 1)]
        self.name_type_func[5][2] = ["GAIN_MODEL_F_N",
                                     f"combo;{self.o.v.GAIN_MODES[params['GAIN_MODEL_F_N']]}",
                                     ";".join(self.o.v.GAIN_MODES), (1, 1)]
        self.name_type_func[6][1] = ["Модель кубсата", "label", "", (1, 1)]
        self.name_type_func[6][2] = ["CUBESAT_MODEL_N",
                                     f"combo;{self.o.v.CUBESAT_MODELS[params['CUBESAT_MODEL_N']]}",
                                     ";".join(self.o.v.CUBESAT_MODELS), (1, 1)]
        self.name_type_func[7][1] = ["Модель чипсата", "label", "", (1, 1)]
        self.name_type_func[7][2] = ["CHIPSAT_MODEL_N",
                                     f"combo;{self.o.v.CHIPSAT_MODELS[params['CHIPSAT_MODEL_N']]}",
                                     ";".join(self.o.v.CHIPSAT_MODELS), (1, 1)]

        # Второй столбец
        self.name_type_func[1][3] = ["Кубсаты", "label", None, (1, 1)]
        self.name_type_func[1][4] = ["CUBESAT_AMOUNT", f"edit;{params['CUBESAT_AMOUNT']}", None, (1, 1)]
        self.name_type_func[2][3] = ["Чипсаты", "label", None, (1, 1)]
        self.name_type_func[2][4] = ["CHIPSAT_AMOUNT", f"edit;{params['CHIPSAT_AMOUNT']}", None, (1, 1)]

        # Третий столбец
        self.name_type_func[1][5] = ["Лобовое сопротивление", "label", None, (1, 1)]
        self.name_type_func[1][6] = ["DYNAMIC_MODEL_aero", f"check;{int(params['DYNAMIC_MODEL_aero'])}", None, (1, 1)]
        self.name_type_func[2][5] = ["Гармоника J₂", "label", None, (1, 1)]
        self.name_type_func[2][6] = ["DYNAMIC_MODEL_j2", f"check;{int(params['DYNAMIC_MODEL_j2'])}", None, (1, 1)]
        self.name_type_func[3][5] = ["Оценка ориентации", "label", None, (1, 1)]
        self.name_type_func[3][6] = ["NAVIGATION_ANGLES", f"check;{int(params['NAVIGATION_ANGLES'])}", None, (1, 1)]
        self.name_type_func[4][5] = ["Разложение сигнала на выходе", "label", None, (1, 1)]
        self.name_type_func[4][6] = ["MULTI_ANTENNA_TAKE", f"check;{int(params['MULTI_ANTENNA_TAKE'])}", None, (1, 1)]
        self.name_type_func[5][5] = ["Разложение сигнала на входе", "label", None, (1, 1)]
        self.name_type_func[5][6] = ["MULTI_ANTENNA_SEND", f"check;{int(params['MULTI_ANTENNA_SEND'])}", None, (1, 1)]
        self.name_type_func[6][5] = ["Навигация (при отладке выключить)", "label", None, (1, 1)]
        self.name_type_func[6][6] = ["IF_NAVIGATION", f"check;{int(params['IF_NAVIGATION'])}", None, (1, 1)]

    def local_save_params(self):
        """Сохранение настроенных параметров"""
        if self.textboxes['Параметры'].text() != "":
            self.o.v.DESCRIPTION = self.textboxes['Параметры'].text()

        self.o.v.dT = float(self.comboboxes['dT'].currentText())
        self.o.v.TIME = float(self.comboboxes['TIME'].currentText())
        self.o.v.START_NAVIGATION = self.comboboxes['START_NAVIGATION_N'].currentText()
        self.o.v.START_NAVIGATION_N = self.o.v.NAVIGATIONS.index(self.o.v.START_NAVIGATION)
        self.o.v.GAIN_MODEL_C = self.comboboxes['GAIN_MODEL_C_N'].currentText()
        self.o.v.GAIN_MODEL_C_N = self.o.v.GAIN_MODES.index(self.o.v.GAIN_MODEL_C)
        self.o.v.GAIN_MODEL_F = self.comboboxes['GAIN_MODEL_F_N'].currentText()
        self.o.v.GAIN_MODEL_F_N = self.o.v.GAIN_MODES.index(self.o.v.GAIN_MODEL_F)
        self.o.v.CUBESAT_MODEL = self.comboboxes['CUBESAT_MODEL_N'].currentText()
        self.o.v.CUBESAT_MODEL_N = self.o.v.CUBESAT_MODELS.index(self.o.v.CUBESAT_MODEL)
        self.o.v.CHIPSAT_MODEL = self.comboboxes['CHIPSAT_MODEL_N'].currentText()
        self.o.v.CHIPSAT_MODEL_N = self.o.v.CHIPSAT_MODELS.index(self.o.v.CHIPSAT_MODEL)

        self.o.v.CUBESAT_AMOUNT = int(self.textboxes['CUBESAT_AMOUNT'].text())
        self.o.v.CHIPSAT_AMOUNT = int(self.textboxes['CHIPSAT_AMOUNT'].text())

        self.o.v.DYNAMIC_MODEL['aero drag'] = self.checkboxes['DYNAMIC_MODEL_aero'].isChecked()
        self.o.v.DYNAMIC_MODEL['j2'] = self.checkboxes['DYNAMIC_MODEL_j2'].isChecked()
        self.o.v.NAVIGATION_ANGLES = self.checkboxes['NAVIGATION_ANGLES'].isChecked()
        self.o.v.MULTI_ANTENNA_TAKE = self.checkboxes['MULTI_ANTENNA_TAKE'].isChecked()
        self.o.v.MULTI_ANTENNA_SEND = self.checkboxes['MULTI_ANTENNA_SEND'].isChecked()
        self.o.v.IF_NAVIGATION = self.checkboxes['IF_NAVIGATION'].isChecked()

        self.o.v.save_params()
        self.config_choose_n = len(self.o.v.config_choose) - 1
        self.close()
        self.initUI()

    def local_load_params(self, i: int):
        self.config_choose_n = i
        self.close()
        self.initUI()

    def local_remove_params(self, i: int):
        if len(self.o.v.config_choose) > 1:
            self.o.v.remove_params(i)
        if i <= self.config_choose_n:
            self.config_choose_n -= 1
        self.close()
        self.initUI()

    def main_run(self):
        """Функция запуска численного моделирования, выводы результатов
        ДОРАБОТАТЬ: какие именно результаты выводить - в чекбокс"""
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
    window = Window(o=o)
    return app, window
