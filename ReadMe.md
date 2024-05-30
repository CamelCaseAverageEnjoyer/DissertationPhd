# Python-проект для аспирантского диплома
| Тема дипломной работы | Определение относительного движения космических аппаратов в групповом полёте |
|-----------------------|------------------------------------------------------------------------------|
| Кафедра               | Математического моделирования и прикладной математики                        |
| Физтех-школа          | Прикладной Математики и Информатики (ФПМИ)                                   |

<img src="source/robot1.png" alt="robot image" width="50"/> "Помнишь, я тебе говорила про мусор, который стоит? Стоит и смердит? Так вот — это была метафора. Я имела в виду тебя."

Плохой пример обнаружения чипсата только по оценке расстояния:

<img src="source/example.gif" width="400">

### Основная часть проекта расположена в каталоге _srs/kiamfemtosat_:
| main       | экспресс-настройка параметров, запуска численного интегрирования и онлайн-навигации |
|------------|------------------------------------------------------------------------------------|
| navigation | решение                                                                            |

- main - файл экспресс-настройки параметров, запуска численного интегрирования и онлайн-навигации
- navigation
- main_objects.py - файл классов:
  - Cosmetic - класс вывода текста (ну и прочих приколюх, которых я пока не придумал)
  - FemtoSat - класс для фемтоспутников
  - CubeSat - класс для чипсатов
  - KalmanFilter - отдельный класс для фильтра Калмана (вообще, можно было бы сделать класс методов навигации, но у меня этих методов пока не очень много)
  - PhysicModel - отдельный класс динамической модели системы
  - Objects - класс, содержащий всё вышеперечисленное
- plot_func.py - файл с функциями отображения
- experiments.py - файл с неосновным численным моделированием. Как правило, здесь функции либо проверяют работоспособность методов, либо строят красивые графики
- tiny_functions.py - файл с микро-функциями из разных областей (округов, резиденций; они развязали великую войну, приведшую к великой скорби)
- frases.py - мне скучно было просто так заниматься дипломом
- helper.py - не обращайте внимания, tmp-файл для быстрых тестовых прогонов чего-либо