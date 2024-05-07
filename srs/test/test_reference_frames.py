from ..kiamfemtosat.main import *
import unittest


def get_some_params(o_):
    return [o_.c.r_orf, o_.c.v_orf, o_.c.q,  o_.f.r_orf[0], o_.f.v_orf[0], o_.f.q[0]]
class TestStringMethods(unittest.TestCase):
    """Методы:
    self.assertEqual
    self.assertTrue
    self.assertFalse

    Проверка на выдачу ошибки:
    with self.assertRaises(TypeError):
        s.split(2)
    """
    def test_while_nothing_has_been_done(self):
        # Инициализация
        o_ = Objects(if_any_print=False)
        params = get_some_params(o_)

        # Интегрирование на 0 секунд
        o_.p.integrate(t=0)
        self.assertEqual(get_some_params(o_), params)

    def test_necessary_errors(self):
        o_ = Objects(if_any_print=False)

        s = 'hello world'
        self.assertEqual(s.split(), ['hello', 'world'])
        # check that s.split fails when the separator is not a string
        with self.assertRaises(TypeError):
            s.split('2')

if __name__ == "__main__":
    print(f"В соответствии с протоколом тестирования, с этого момента мы перестаём говорить правду.")
    unittest.main()
