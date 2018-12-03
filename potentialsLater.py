import numpy as np
import cmd
import config as c

'''
Inspired from https://docs.python.org/3.3/library/cmd.html
'''


def test():
    print("It works")


# Defining variables
omega = np.sqrt(c.ko / c.mass)


class PotClass(cmd.Cmd):
    intro = \
        'This program will use the Crank-Nicolson to propegate a Gaussian '\
        'wavepacket in \n various potential functions. Type help or ? for a' \
        ' of commands. \n'
    prompt = '(CN-Gauss) '
    file = None

    def do_free(self, arg):
        'A zero potential function. Free particle'
        def pot():
            val = 0
            return val

    def do_SHO(self, arg):
        'The potential energy function for the simple harmonic oscillator.'
        def pot(x):
            val = (1 / 2) * c.mass * (omega ** 2) * (x ** 2)
            return val

    def do_MM(self, arg):
        'The potential used in the adv. math methods problem.'
        def pot(Vo, x1, a, x):
            val = Vo * np.exp(
                -((x - x1) ** 2)
                / 2 * (a ** 2)
            )
            return val

    def do_quit(self, arg):
        'Close the program.'
        print('Thank you!')
        self.close()
        return True

    def close(self):
        if self.file:
            self.file.close()
            self.file = None


# if __name__ == '__main__':
#    PotClass().cmdloop()

potfuncs = PotClass()
