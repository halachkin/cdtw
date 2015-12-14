import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from cydtw import *


class dtw(cydtw):

    def __init__(self, ref, query, settings=Settings()):
        cydtw.__init__(self, ref, query, settings)
        self._ref = ref
        self._query = query

    def plot_mat_path(self, colormap=cm.Greys_r):
        ils = []
        jls = []
        for i, j in self.get_path():
            ils.append(i), jls.append(j)

        fig = plt.figure(1)

        ax = fig.add_subplot(111)
        ax.plot(jls, ils, 'w', lw=1.5)
        cax = ax.imshow(self._cost, interpolation='nearest',
                        cmap=colormap,
                        origin='lower')

        ax.set_xlim([-0.5, self._cost.shape[1] - 0.5])
        ax.set_ylim([self._cost.shape[0] - 0.5, -0.5])

        # divider = make_axes_locatable(cax)
        # cax = divider.append_axes("right", size="5%", pad=0.2)

        cbar = plt.colorbar(cax)
        # plt.show()

        return ax

    def plot_seq_mat_path(self, colormap=cm.Greys_r):

        plt.figure(1)
        pr = 6

        ax1 = plt.subplot2grid((pr, pr), (0, 1), colspan=4)  # query
        ax2 = plt.subplot2grid((pr, pr), (1, 0), rowspan=4)  # ref
        ax3 = plt.subplot2grid((pr, pr), (1, 1), colspan=4,
                               rowspan=4,
                               sharex=ax1,
                               sharey=ax2)  # cost and path

        # axis limit
        # axis((xmin,xmax,ymin,ymax))
        ax1.axis((0, self._cost.shape[1], min(self._query), max(self._query)))

        ax2.axis((np.min(self._ref), np.max(self._ref),
                  0, self._cost.shape[0]))

        ax3.axis((-0.5, self._cost.shape[1] - 0.5,
                  self._cost.shape[0] - 0.5, -0.5))

        ax1.plot(self._query, 'k')
        ax2.plot(self._ref, np.arange(len(self._ref)), 'k')

        ax3.imshow(self._cost,
                   interpolation='nearest',
                   cmap=colormap,
                   aspect="auto",
                   origin='lower')

        ils = []
        jls = []
        for i, j in self.get_path():
            ils.append(i), jls.append(j)

        ax3.plot(jls, ils, 'w', lw=1.5)

        ax3.xaxis.set_visible(False)
        ax2.yaxis.set_visible(False)
        ax2.xaxis.set_visible(False)
        ax1.yaxis.set_visible(False)

        # for i, row in enumerate(self._cost):
        #     for j, c in enumerate(row):
        #         if c<10:
        #             ax3.text(j-.1, i+.1, int(c), fontsize=13)
        #         else:
        #             ax3.text(j-.2, i+.1, int(c), fontsize=13)

        return ax1, ax2, ax3

    def plot_alignment(self, coef=2, title=''):
        copy_ref = np.copy(self._ref)

        for i in range(len(copy_ref)):
            copy_ref[i] += coef * np.max(self._ref)

        plt.figure(1)
        plt.plot(copy_ref, lw=3)
        plt.plot(self._query, lw=3)
        plt.xlim(-0.5, max(len(self._ref), len(self._query)) + 1)
        for val in self._path:
            plt.plot([val[0], val[1]], \
                     [copy_ref[val[0]], self._query[val[1]]], 'k', lw=0.5)

        plt.axis('off')
        plt.show()


if __name__ == '__main__':
    d = dtw(np.random.rand(10), np.random.rand(
        10), Settings(compute_path=True))

    d.plot_seq_mat_path()
    plt.show()
