from __future__ import absolute_import
from RateDistribution import Distro
from Queue import Empty
import pickle
from numpy.random import get_state, set_state, seed
from unittest import TestCase


class TestDistro(TestCase):
    temp_pickle_file = 'pickles/temp_distro.p'


    def setUp(self):
        seed()
        self.seed = get_state()
        self.distro = Distro('exp', size=3, scale=0.5)

    def tearDown(self):
        import os
        try:
            os.remove(temp_pickle_file)
        except OSError:
            pass
        # not sure if this works...
        del self.distro


    def test_distro_dump_persistence(self):
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
                'Dumping the distro into a pickle altered the queue content'
                )

    def test_distro_load_equality(self):
        """
        Makes sense to test only after test_distro_dump_persistence
        succeeded.
        """
        # dump distro
        # load it as distro2
        # compare queue size and content
        # compare draw_fct
        # compare even __dict__ ?
