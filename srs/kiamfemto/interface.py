import sys
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from my_plot import *
from simulation import save_simulation_trajectories, load_simulation_trajectories

ICON_SIZE = 70

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

        # Подгрузка параметров + инициализация
        self.o.v.load_params(i=self.config_choose_n)
        self.o.init_classes()

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
                if ";" in t:
                    s = int(t.split(";")[1])
                    label.setPixmap(pix.scaled(s, s, Qt.AspectRatioMode.KeepAspectRatio,
                                               Qt.TransformationMode.FastTransformation))
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
        y_all = 0
        n = 0

        # >>>>>>>>>>>> Кнопки <<<<<<<<<<<<
        y = 0
        self.name_type_func[y][n] = [self.path + "robot1.png", "button", talk, (1, 1)]  # Поболтать
        y += 1
        self.name_type_func[y][n] = [self.path + "integral.png", "button", self.main_run, (1, 1)]  # Моделирование
        y += 1
        self.name_type_func[y][n] = [self.path + "antenna.png", "button", plot_model_gain, (1, 1)]  # Диаграммы
        y += 1
        self.name_type_func[y][n] = [self.path + "air.png", "button", plot_atmosphere_models, (1, 1)]
        y += 1
        self.name_type_func[y][n] = [self.path + "animation.png", "button", animate_reference_frames, (1, 1)]
        y += 1
        self.name_type_func[y][n] = [self.path + "plot.png", "button", lambda x=self.o: plot_distance(x), (1, 1)]
        y += 1
        self.name_type_func[y][n] = [self.path + "param.png", "button", self.plot_1_param, (1, 1)]
        y += 1
        self.name_type_func[y][n] = [self.path + "orbit.png", "button", lambda x=self.o: plot_all(x), (1, 1)]
        y += 1
        self.name_type_func[y][n] = [self.path + "signal.png", "button", lambda x=self.o: plot_signals(x), (1, 1)]
        y += 1
        self.name_type_func[y][n] = [self.path + "save.png", "button", self.local_save_trajectories, (1, 1)]
        y += 1
        self.name_type_func[y][n] = [self.path + "load.png", "button", self.local_load_trajectories, (1, 1)]
        y += 1
        self.name_type_func[y][n] = [self.path + "eraser.png", "button", self.local_remove_trajectories, (1, 1)]
        y += 1
        n += 1

        # >>>>>>>>>>>> Параметры численного моделирования <<<<<<<<<<<<

        # Первый столбец
        y = 1
        self.name_type_func[y][n+0] = ["Шаг по времени dT", "label", None, (1, 1)]
        self.name_type_func[y][n+1] = ["dT", f"combo;{params['dT']}", ";".join(self.o.v.dTs), (1, 1)]
        y += 1
        self.name_type_func[y][n+0] = ["Время интегрирования T", "label", "", (1, 1)]
        self.name_type_func[y][n+1] = ["TIME", f"combo;{params['TIME']}", ";".join(self.o.v.Ts), (1, 1)]
        y += 1
        self.name_type_func[y][n+0] = ["Априорная информация", "label", "", (1, 1)]
        self.name_type_func[y][n+1] = ["START_NAVIGATION_N",
                                       f"combo;{self.o.v.NAVIGATIONS[params['START_NAVIGATION_N']]}",
                                       ";".join(self.o.v.NAVIGATIONS), (1, 1)]
        y += 1
        self.name_type_func[y][n+0] = ["Антенны кубсатов", "label", "", (1, 1)]
        self.name_type_func[y][n+1] = ["GAIN_MODEL_C_N",
                                       f"combo;{self.o.v.GAIN_MODES[params['GAIN_MODEL_C_N']]}",
                                       ";".join(self.o.v.GAIN_MODES), (1, 1)]
        y += 1
        self.name_type_func[y][n+0] = ["Антенны чипсатов", "label", "", (1, 1)]
        self.name_type_func[y][n+1] = ["GAIN_MODEL_F_N",
                                       f"combo;{self.o.v.GAIN_MODES[params['GAIN_MODEL_F_N']]}",
                                       ";".join(self.o.v.GAIN_MODES), (1, 1)]
        y += 1
        self.name_type_func[y][n+0] = ["Модель кубсата", "label", "", (1, 1)]
        self.name_type_func[y][n+1] = ["CUBESAT_MODEL_N",
                                       f"combo;{self.o.v.CUBESAT_MODELS[params['CUBESAT_MODEL_N']]}",
                                       ";".join(self.o.v.CUBESAT_MODELS), (1, 1)]
        y += 1
        self.name_type_func[y][n+0] = ["Модель чипсата", "label", "", (1, 1)]
        self.name_type_func[y][n+1] = ["CHIPSAT_MODEL_N",
                                       f"combo;{self.o.v.CHIPSAT_MODELS[params['CHIPSAT_MODEL_N']]}",
                                       ";".join(self.o.v.CHIPSAT_MODELS), (1, 1)]
        y += 1
        n += 2
        y_all = max(y_all, y)


        # Второй столбец
        y = 1
        self.name_type_func[y][n+0] = ["Кубсаты", "label", None, (1, 1)]
        self.name_type_func[y][n+1] = [self.path + f"{self.o.v.CUBESAT_MODELS[params['CUBESAT_MODEL_N']]}.png",
                                       f"image;{ICON_SIZE}", None, (1, 1)]
        self.name_type_func[y][n+2] = ["CUBESAT_AMOUNT", f"edit;{params['CUBESAT_AMOUNT']}", None, (1, 1)]
        y += 1
        self.name_type_func[y][n+0] = ["Чипсаты", "label", None, (1, 1)]
        self.name_type_func[y][n+1] = [self.path + f"chipsat.png", f"image;{ICON_SIZE}", None, (1, 1)]
        self.name_type_func[y][n+2] = ["CHIPSAT_AMOUNT", f"edit;{params['CHIPSAT_AMOUNT']}", None, (1, 1)]
        y += 1
        self.name_type_func[y][n+0] = ["q", "label", None, (1, 1)]
        self.name_type_func[y][n+1] = ["q", f"edit;{params['q']}", None, (1, 2)]
        y += 1
        self.name_type_func[y][n+0] = ["p", "label", None, (1, 1)]
        self.name_type_func[y][n+1] = ["p", f"edit;{params['p']}", None, (1, 2)]
        y += 1
        self.name_type_func[y][n+0] = ["r", "label", None, (1, 1)]
        self.name_type_func[y][n+1] = ["r", f"edit;{params['r']}", None, (1, 2)]
        y += 1
        n += 3
        y_all = max(y_all, y)

        # Третий столбец
        y = 1
        self.name_type_func[y][n+0] = ["Лобовое сопротивление", "label", None, (1, 1)]
        self.name_type_func[y][n+1] = ["DYNAMIC_MODEL_aero", f"check;{int(params['DYNAMIC_MODEL_aero'])}", None, (1, 1)]
        y += 1
        self.name_type_func[y][n+0] = ["Гармоника J₂", "label", None, (1, 1)]
        self.name_type_func[y][n+1] = ["DYNAMIC_MODEL_j2", f"check;{int(params['DYNAMIC_MODEL_j2'])}", None, (1, 1)]
        y += 1
        self.name_type_func[y][n+0] = ["Оценка ориентации", "label", None, (1, 1)]
        self.name_type_func[y][n+1] = ["NAVIGATION_ANGLES", f"check;{int(params['NAVIGATION_ANGLES'])}", None, (1, 1)]
        y += 1
        self.name_type_func[y][n+0] = ["Разложение сигнала на выходе", "label", None, (1, 1)]
        self.name_type_func[y][n+1] = ["MULTI_ANTENNA_TAKE", f"check;{int(params['MULTI_ANTENNA_TAKE'])}", None, (1, 1)]
        y += 1
        self.name_type_func[y][n+0] = ["Разложение сигнала на входе", "label", None, (1, 1)]
        self.name_type_func[y][n+1] = ["MULTI_ANTENNA_SEND", f"check;{int(params['MULTI_ANTENNA_SEND'])}", None, (1, 1)]
        y += 1
        self.name_type_func[y][n+0] = ["Навигация (при отладке выключить)", "label", None, (1, 1)]
        self.name_type_func[y][n+1] = ["IF_NAVIGATION", f"check;{int(params['IF_NAVIGATION'])}", None, (1, 1)]
        y += 1
        n += 2
        y_all = max(y_all, y)

        self.name_type_func[0][1] = ["Параметры численного моделирования", "label", None, (1, n - 1)]

        # >>>>>>>>>>>> Автоматическое заполнение сохранений <<<<<<<<<<<<
        y_save = len(self.o.v.config_choose)
        self.name_type_func[y_save + 1][n+0] = ["Параметры", "edit;Название", None, (1, 2)]
        self.name_type_func[y_save + 1][n+2] = ["Сохранить", "button", self.local_save_params, (1, 1)]
        self.name_type_func[y_save + 2][n+0] = ["Применить введённые параметры", "button", self.apply_params, (1, 3)]

        for i in range(y_save):
            self.name_type_func[i+1][n+0] = [data.iloc[i, 0], "label", None, (1, 1)]
            self.name_type_func[i+1][n+1] = [f"Загрузить", "button", lambda j=i: self.local_load_params(i=j),
                                                   (1, 1)]
            self.name_type_func[i+1][n+2] = [f"Удалить", "button", lambda j=i: self.local_remove_params(i=j),
                                                   (1, 1)]

        # >>>>>>>>>>>> Отрисовка <<<<<<<<<<<<
        # Диаграмма разложения сигналов
        n = 1
        self.name_type_func[y_all+0][n+0] = [f"Разложения сигналов: {params['MULTI_ANTENNA_SEND']} -> "
                                              f"{params['MULTI_ANTENNA_TAKE']}", "label", None, (1, 3)]
        self.name_type_func[y_all+1][n+0] = [self.path + "send0_take0.png", "image", None, (3, 3)]

        # Текст о долготе полёта
        self.name_type_func[y_all][n+3] = [self.o.time_message(params['TIME']), "label", None, (1, 3)]

    def apply_params(self):
        """Применение настроенных параметров. Должно быть согласовано с config.get_saving_params"""
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
        self.o.v.KALMAN_COEF['q'] = [float(self.textboxes['q'].text())] * 2
        self.o.v.KALMAN_COEF['p'] = [float(self.textboxes['p'].text())] * 4
        self.o.v.KALMAN_COEF['r'] = float(self.textboxes['r'].text())

        self.o.v.DYNAMIC_MODEL['aero drag'] = self.checkboxes['DYNAMIC_MODEL_aero'].isChecked()
        self.o.v.DYNAMIC_MODEL['j2'] = self.checkboxes['DYNAMIC_MODEL_j2'].isChecked()
        self.o.v.NAVIGATION_ANGLES = self.checkboxes['NAVIGATION_ANGLES'].isChecked()
        self.o.v.MULTI_ANTENNA_TAKE = self.checkboxes['MULTI_ANTENNA_TAKE'].isChecked()
        self.o.v.MULTI_ANTENNA_SEND = self.checkboxes['MULTI_ANTENNA_SEND'].isChecked()
        self.o.v.IF_NAVIGATION = self.checkboxes['IF_NAVIGATION'].isChecked()

        my_print('Параметры применены!', color='c')

    def local_save_params(self):
        """Сохранение настроенных параметров"""
        self.apply_params()
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

    def local_save_trajectories(self):
        text, ok = QInputDialog.getText(self, 'Сохранение траектории', 'Введите название файла:')
        if ok:
            save_simulation_trajectories(o=self.o, text=f"{self.o.v.path_sources}trajectories/{text}")

    def local_load_trajectories(self):
        path = f"{self.o.v.path_sources}trajectories/"
        text = QFileDialog.getOpenFileName(self, 'Open file', path)[0]
        if text != "":
            load_simulation_trajectories(o=self.o, text=path + text.split(path)[1])

    def local_remove_trajectories(self):
        from os import listdir, remove
        items = sorted([s for s in listdir(f"{self.o.v.path_sources}trajectories")])
        if len(items) > 0:
            text, ok = QInputDialog.getItem(self, "Удаление", "Выберите траекторию для удаления", items, 0, False)
            if ok:
                remove(f"{self.o.v.path_sources}trajectories/{text}")

    def plot_1_param(self):
        d = self.o.p.record
        items = sorted(d.columns)

        # ручная сортировка
        items = [i for i in items if not ('orf' in i or 'irf' in i)] + \
                [i for i in items if 'orf' in i] + [i for i in items if 'irf' in i]

        y_s = []
        ok = True
        while ok:
            text, ok = QInputDialog.getItem(self, "Отображение", "Выберите насколько параметров", items, 0, False)
            y_s.append(text)
        for y in y_s[:-1]:  # Костыль - в конце добавляется items[0]
            plt.plot(d['t'].to_list(), d[y].to_list(), label=y)
        plt.legend()
        plt.grid()
        plt.show()

    def main_run(self):
        """Функция запуска численного моделирования, выводы результатов"""
        self.o.integrate(t=self.o.v.TIME, animate=False)

        # Вывод результатов
        tmp = np.array(self.o.p.record[f'{self.o.f.name} KalmanPosError r {0}'].to_list())  # Для чипсата id=0
        print(f"Математическое ожидание ошибки: {tmp.mean():.2f} м, Среднее отклонение ошибки: {tmp.std():.2f} м")
        talk_decision(cnd=self.o.v.IF_TALK)
        plot_distance(self.o)

def interface_window(o):
    app = QApplication(sys.argv)
    window = Window(o=o)
    return app, window
