import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

from cydtw import cydtw, Settings


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


        # To display two dimensional array as one dimensional string, we need to convert 2-dim points to 1 dims points
        # by taking the norm of each points

        if len(self._query.shape) == 1:
            query_1d = self._query
            ref_1d = self._ref
        else:
            query_1d = np.sum(np.abs(self._query)**2,axis=-1)**(1./2)
            ref_1d = np.sum(np.abs(self._ref) ** 2, axis=-1) ** (1. / 2)

        ax1.axis((0, len(query_1d), min(query_1d), max(query_1d)))
        ax2.axis((np.min(ref_1d), np.max(ref_1d), 0, len(ref_1d)))

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

    # Two dimensional
    arr1 = np.array([[1, 0, 0], [5, 0, 0], [4, 0, 0], [2, 0, 0]])
    arr2 = np.array([[1, 0, 0], [2, 0, 0], [4, 0, 0], [1, 0, 0], [2, 0, 0]])
    settings = Settings(compute_path=True)
    d = dtw(arr1, arr2, settings)
    print(d.get_cost())
    print(d.get_path())
    d.plot_seq_mat_path()
    plt.show()

    # One dimensional
    arr1 = np.array([1, 5, 4, 2])
    arr2 = np.array([1, 2, 4, 1])
    settings = Settings(compute_path=True)
    d = dtw(arr1, arr2, settings)
    print(d.get_cost())
    print(d.get_path())
    d.plot_seq_mat_path()
    plt.show()
