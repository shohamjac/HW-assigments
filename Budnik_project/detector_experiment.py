import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


class DetectorExperiment:
    """
    A class for scope acquisition with many excitations.
    """

    def __init__(self, data: pd.DataFrame, number_of_detections=1000, pts_per_trig=1500, trig_point=700):
        self.number_of_detections = number_of_detections
        self.pts_per_trig = pts_per_trig
        self.data = data
        self.dt = data.second[1] - data.second[0]
        self.energy_list = None
        self.steepness_list = None
        self.trig_point = trig_point

    def get_energy_list(self, **kwargs):
        # if self.energy_list is None:
        #     self.energy_list = [self.get_excitation_energy_with_edge_detection(excitation, **kwargs) for excitation in
        #                         range(self.number_of_detections)]
        # return self.energy_list
        return [self.get_excitation_energy_with_edge_detection(excitation, **kwargs) for excitation in
                                range(self.number_of_detections)]

    def get_steepness_list(self, length, **kwargs):
        # if self.steepness_list is None:
        #     self.steepness_list = [self.get_excitation_steepness(excitation, length, **kwargs) for excitation in
        #                            range(self.number_of_detections)]
        # return self.steepness_list
        return [self.get_excitation_steepness(excitation, length, **kwargs)
                for excitation in range(self.number_of_detections)]

    def get_width_list_based_on_peak(self, peak_multiplier, **kwargs):
        edge_val = [peak_multiplier * min(self.data.Volt[excitation * self.pts_per_trig:
                                                         (excitation + 1) * self.pts_per_trig])
                    for excitation in range(self.number_of_detections)]
        startandstop = [self.get_excitation_start_and_stop(excitation, edge_val=edge_val[excitation], **kwargs)
                        for excitation in range(self.number_of_detections)]
        return [self.dt*(tup[1] - tup[0]) for tup in startandstop]

    def get_time_list(self, **kwargs):
        # if self.steepness_list is None:
        #     self.steepness_list = [self.get_excitation_steepness(excitation, length, **kwargs) for excitation in
        #                            range(self.number_of_detections)]
        # return self.steepness_list
        return [self.get_excitation_time(excitation, **kwargs)
                for excitation in range(self.number_of_detections)]

    def get_excitation_energy(self, exc_number, start, stop):
        """
        Integrating the excitation signal from 'start' to 'end'
        :param exc_number: the excitation number in the data
        :param start: the start of the integration point
        :param stop: the stop of the integration point
        :return: The integration value in unit of Volt*Sec
        """
        return self.dt * sum(
            self.data.Volt[exc_number * self.pts_per_trig + start:exc_number * self.pts_per_trig + stop])

    def get_excitation_energy_with_edge_detection(self, exc_number, **kwargs):
        """
        Returns the energy from the start point to the last point where the signal is smaller than end_edge
        :param RMS_multiplier: In the case edge_val is None, we set it to be RMS*RMS_multiplier.
        :param exc_number: the excitation number in the data
        :param edge_val: the value for the end detection
        :param start: The start of the integration point
        :param plot: produce a plot with the end and start point. mainly for debug reasons.
        :return: The integration value in unit of Volt*Sec
        """
        start, stop = self.get_excitation_start_and_stop(exc_number, **kwargs)
        return self.dt * sum(
            self.data.Volt[exc_number * self.pts_per_trig + start:exc_number * self.pts_per_trig + stop])

    def get_excitation_time(self, exc_number):
        return self.data.second[exc_number * self.pts_per_trig + self.trig_point]

    def get_excitation_steepness(self, exc_number, length, **kwargs):
        #print(exc_number)
        start, _ = self.get_excitation_start_and_stop(exc_number, **kwargs)
        return self.dt * sum(
            self.data.Volt[exc_number * self.pts_per_trig + start:exc_number * self.pts_per_trig + start + int(length/self.dt)])

    def get_excitation_start_and_stop(self, exc_number, edge_val=None, plot=False, RMS_multiplier=3):
        '''

        :param exc_number:
        :param edge_val: negative value in which we cut by
        :param plot:
        :param RMS_multiplier:
        :return:
        '''
        if edge_val is None:
            edge_val = - RMS_multiplier * np.std(self.data.Volt[exc_number * self.pts_per_trig
                                                                :exc_number * self.pts_per_trig + int(
                self.trig_point / 2)])
        min_voltage_ind = np.argmin(self.data.Volt[exc_number * self.pts_per_trig:(exc_number + 1) * self.pts_per_trig])
        start = np.asarray(np.array(self.data.Volt[exc_number * self.pts_per_trig
                                                   :exc_number * self.pts_per_trig + min_voltage_ind])
                           > edge_val).nonzero()[0][-1]
        try:  # It’s easier to ask for forgiveness than permission
            stop = np.asarray(np.array(self.data.Volt[exc_number * self.pts_per_trig
                                                      + min_voltage_ind:(exc_number + 1) * self.pts_per_trig])
                              > edge_val).nonzero()[0][0] + min_voltage_ind
        except IndexError:
            stop = self.pts_per_trig - 1
        if plot:
            plt.plot(self.data.second[exc_number * self.pts_per_trig:(exc_number + 1) * self.pts_per_trig],
                     self.data.Volt[exc_number * self.pts_per_trig:(exc_number + 1) * self.pts_per_trig], 'b-')
            plt.axvline(self.data.second[exc_number * self.pts_per_trig + start], color='r')
            #plt.axvline(self.data.second[exc_number * self.pts_per_trig + start + 40], color='g')
            plt.axvline(self.data.second[exc_number * self.pts_per_trig + stop], color='r')
            plt.axhline(edge_val, color='r')

            plt.show()
        return start, stop
