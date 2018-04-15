from __future__ import absolute_import
from ..RateDistribution import Distro
from Queue import Empty
import pickle
from numpy.random import get_state, set_state, seed
from unittest import TestCase


class TestDistro(TestCase):
    temp_pickle_file = 'temp_distro.p'


    def setUp(self):
        seed()
        self.queue_size = 3
        self.seed = get_state()
        self.distro = Distro('exp', size=self.queue_size, scale=0.5)

    def tearDown(self):
        import os
        try:
            os.remove(self.temp_pickle_file)
        except OSError:
            pass
        # not sure if this works...
        del self.distro


    def test_distro_queue_dump(self):
        """
        Integrity of Distro.queue when exporting the Distro object
        to a pickle.
        """
        import pickle
        queue_list = []
        queue_size = self.distro.queue.qsize()
        while True:
            try:
                queue_list.append(self.distro.queue.get_nowait())
            except Empty:
                break
        for el in queue_list:
            self.distro.queue.put_nowait(el)
        with open(self.temp_pickle_file, 'wb') as f:
            pickle.dump(self.distro, f)
        queue_list_after = []
        queue_size_after = self.distro.queue.qsize()
        while True:
            try:
                queue_list_after.append(self.distro.queue.get_nowait())
            except Empty:
                break
        self.assertEqual(
                queue_size,
                queue_size_after,
                'Dumping the distro into a pickle altered the queue length'
                )
        self.assertEqual(
                queue_list,
                queue_list_after,
                'Dumping the distro into a pickle altered the queue' \
                    'content:\n{0}'.format('\n'.join(
                        map(str, zip(queue_list, queue_list_after))
                        )
                    )
                )

    def test_distro_queue_load(self):
        """
        Integrity of Distro.queue when loading a Distro object
        form a pickle.

        Makes sense to test only after test_distro_dump_persistence
        succeeded.
        """
        import pickle
        queue_list = []
        queue_size = self.distro.queue.qsize()
        with open(self.temp_pickle_file, 'wb') as f:
            pickle.dump(self.distro, f)
        with open(self.temp_pickle_file, 'rb') as f:
            distro2 = pickle.load(f)
        queue_list= []
        queue_size= self.distro.queue.qsize()
        while True:
            try:
                queue_list.append(self.distro.queue.get_nowait())
            except Empty:
                break
        queue_list_loaded = []
        queue_size_loaded = distro2.queue.qsize()
        while True:
            try:
                queue_list_loaded.append(distro2.queue.get_nowait())
            except Empty:
                break
        self.assertEqual(
                queue_size,
                queue_size_loaded,
                'Dumping the distro into a pickle altered the queue length'
                )
        self.assertEqual(
                queue_list,
                queue_list_loaded,
                'Dumping the distro into a pickle altered the queue' \
                    'content:\n{0}'.format('\n'.join(
                        map(str, zip(queue_list, queue_list_loaded))
                        )
                    )
                )

    def test_distro_v_get_integrity(self):
        """
        Distro.v_get equivalent to repeated Distro.get_val
        """
        pass

    def test_distro_draw_fct_load(self):
        """
        Reproducibility of Distro.draw_fct upon pickleing
        """
        from numpy import array as narray
        # save distro
        with open(self.temp_pickle_file, 'wb') as f:
            pickle.dump(self.distro, f)
        # reset seed with self.seed
        set_state(self.seed)
        # draw 2*self.queue_size>distro_values
        _shape_object = [0]*(2*self.queue_size)
        distro_values = self.distro.v_get(_shape_object)
        # load as distro2
        with open(self.temp_pickle_file, 'rb') as f:
            distro2 = pickle.load(f)
        # reset seed with self.seed
        set_state(self.seed)
        # draw 2*self.queue_size form distro2>distro2_values
        distro2_values = distro2.v_get(_shape_object)
        # assertEqual distro_values and distro2_values
        self.assertEqual(
                distro_values,
                distro2_values,
                'Loading Distro from pickle alters Distro.draw_fct'
                )

