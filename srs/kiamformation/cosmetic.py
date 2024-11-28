"""Переделать файл перед отправкой на РИД"""

def my_print(txt: any, color: str = None, if_print: bool = True, bold: bool = False, if_return: bool = False,
             end: str = '\n') -> None:
    import colorama
    """Функция вывода цветного текста
    :param txt: Выводимый текст
    :param color: Цвет текста {b, g, y, r, c, m}
    :param if_print: Флаг вывода для экономии места
    :param bold: Жирный текст
    :param if_return: Надо ли возвращать строку"""
    color_bar = {"b": colorama.Fore.BLUE, "g": colorama.Fore.GREEN, "y": colorama.Fore.YELLOW, "r": colorama.Fore.RED,
                 "c": colorama.Fore.CYAN, "m": colorama.Fore.MAGENTA, None: colorama.Style.RESET_ALL}
    _txt = f"\033[1m{txt}\033[0m" if bold else txt
    anw = color_bar[color] + f"{_txt}" + colorama.Style.RESET_ALL
    if if_print and color in color_bar.keys():
        print(anw, end=end)
    if if_return:
        return anw

def real_workload_time(n: int, n_total: int, time_begin, time_now) -> str:
    n_remain = n_total - n
    return f"время: {time_now - time_begin}, оставшееся время: {(time_now - time_begin) * n_remain / n}"

def rand_txt() -> str:
    from random import choice
    with open("data/phrases.txt") as f:
        lines = f.readlines()
    return choice(lines).strip()

def talk_aloud(txt):
    from os import remove
    from playsound3 import playsound
    from gtts import gTTS

    s = gTTS(txt, lang='ru')
    s.save('talk_file.mp3')
    playsound('talk_file.mp3')
    remove('talk_file.mp3')

def talk(aloud=True):
    txt = rand_txt()
    my_print(txt, color='c')
    if aloud:
        talk_aloud(txt)

def ending(n: int) -> str:
    if (n % 10) == 1:
        return " "
    if ((n % 10) > 1) and ((n % 10) < 5):
        return "а"
    if ((n % 10) >= 5) or ((n > 9) and (n < 21)):
        return "ов"
    return "ов"
