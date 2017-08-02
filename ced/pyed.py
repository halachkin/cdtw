import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np

from cyed import cyed, Settings


class EditDistance(cyed):

    def __init__(self, ref, query, args, settings=Settings()):
        cyed.__init__(self, ref, query, args, settings)
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


class Dtw(EditDistance):
    def __init__(self, ref, query, settings=Settings()):
        settings.qtse.set_type('no_quantisation')
        EditDistance.__init__(self, ref, query, {}, settings)


class Edr(EditDistance):
    def __init__(self, ref, query, args, settings=Settings()):
        settings.qtse.set_type('edr')
        settings.step.set_type('dp2edr')
        cyed.__init__(self, ref, query, args, settings)


class Erp(EditDistance):
    def __init__(self, ref, query, args, settings=Settings()):
        settings.qtse.set_type('no_quantisation')
        settings.step.set_type('dp2erp')
        cyed.__init__(self, ref, query, args, settings)


if __name__ == '__main__':

    # Two dimensional
    # arr1 = np.array([[4, 0, 0], [4, 0, 0], [2, 0, 0], [4, 0, 0]])
    # arr2 = np.array([[5, 0, 0], [4, 0, 0], [5, 0, 0], [6, 0, 0], [4, 0, 0]])
    # tolerance = np.array([1, 0, 0])
    # settings = Settings(compute_path=True, dist='edr')
    # d = ed(arr1, arr2, tolerance, settings)
    # print(d.get_cost())
    # print(d.get_path())

    # r = np.array([[4, 0, 0], [4, 0, 0], [2, 0, 0], [4, 0, 0]])
    # q = np.array([[5, 0, 0], [4, 0, 0], [5, 0, 0], [6, 0, 0], [4, 0, 0]])
    # tolerance = 1

    r = np.array([4, 4, 2, 4])
    q = np.array([5, 4, 5, 6, 4])
    tolerance = np.array([1])
    gap = np.array([0])

    args = {'sigmas': tolerance, 'gap': gap}

    d = Erp(r, q, args, settings=Settings(
                              dist='euclid',
                              # step='p0sym',  # Sakoe-Chiba symmetric step with slope constraint p = 0
                              window='palival',  # type of the window
                              param=2.0,  # window parameter
                              norm=False,  # normalization
                              compute_path=True))
    print(d.get_cost())
    print(d.get_path())

    # d.plot_seq_mat_path()
    # plt.show()

    # # One dimensional
    # arr1 = np.array([4, 4, 2, 4])
    # arr2 = np.array([5, 4, 5, 6, 4])
    # tolerance = np.array([1])
    # settings = Settings(compute_path=True)
    # d = ed(arr1, arr2, tolerance, settings)
    # print(d.get_cost())
    # print(d.get_path())
    # d.plot_seq_mat_path()
    # plt.show()
