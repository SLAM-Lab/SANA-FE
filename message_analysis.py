import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import math
import os
import enum
import logging
import sys

np.set_printoptions(threshold=100000)
np.seterr(invalid='raise')

FORMAT = "[%(funcName)s:%(lineno)d] %(message)s"
from functools import partial, partialmethod
logging.TRACE = 5
logging.addLevelName(logging.TRACE, "TRACE")
logging.Logger.trace = partialmethod(logging.Logger.log, logging.TRACE)
logging.trace = partial(logging.log, logging.TRACE)
#logging.basicConfig(format=FORMAT, level=logging.TRACE, stream=sys.stdout)
#logging.basicConfig(format=FORMAT, level=logging.INFO, stream=sys.stdout)
logging.basicConfig(format=FORMAT, level=logging.WARNING, stream=sys.stdout)


def hw_str_to_core(hw_str):
    tile, core = hw_str.split('.')
    core = int(core)
    tile = int(tile)
    return tile, core


def graph_index(x, y, link):
    assert(x < 8)
    assert(y < 4)
    assert(link < 12)
    return (x*4*12) + (y*12) + link


def reverse_graph_index(index):
    x = index // (4*12)
    y = (index % (4*12)) // 12
    link = index % 12
    assert(x < 8)
    assert(y < 4)
    assert(link < 12)
    return (x, y, link)


def calculate_a_k(k, arrival_rate, service_distribution_pdf):
    x = service_distribution_pdf[0]
    y = ((arrival_rate * x) ** k) / math.factorial(k)
    y = y * np.exp(-arrival_rate * x) * service_distribution_pdf[1]

    #print(f"arrival:{arrival_rate} k:{k} dist:{service_distribution_pdf[1]}")
    #print(f"x:{x} y:{y}")

    a = np.sum(y)
    #print(f"Probability of k messages (a_k):{a}")
    # a is a probability and must be < 1
    assert(a <= 1.0)

    # Probability that exactly k messages arrive during a service time
    return a


def calculate_pi(K, arrival_rate, service_distribution_pdf):
    pi_prime = np.zeros((K))
    pi = np.zeros((K))

    # Calculate the arrival (a_k) probabilities
    a = np.zeros((K))
    for k in range(0, K):
        a[k] = calculate_a_k(k, arrival_rate, service_distribution_pdf)

    #print(f"Probability of k messages arriving: {a}")
    pi_prime[0] = 1
    for k in range(0, K-1):
        # Compute every value of pi'
        s = 0
        for j in range(1, k+1):
            s += pi_prime[j] * a[k-j+1]
        pi_prime[k+1] = (1.0/a[0]) * (pi_prime[k] - s - a[k])
    #print(f"pi_prime:{pi_prime}")

    pi[0] = 1.0 / np.sum(pi_prime)
    for k in range(1, K):
        pi[k] = pi[0] * pi_prime[k]

    #print(f"pi[0]:{pi}")
    return pi


def create_pdf(samples):
    if len(samples) == 0:
        return None
        #return np.array(0), np.array(0)

    samples.sort()
    counts = []
    values = []
    prev = None
    for v in samples:
        if v == prev:
            counts[-1] += 1
        else:
            values.append(v)
            counts.append(1)
        prev = v

    values = np.array(values)
    counts = np.array(counts)
    probabilities = np.array(counts) / np.sum(counts)

    # Total probability should roughly=1
    assert abs(np.sum(probabilities) - 1.0) < 0.001
    assert values.shape == probabilities.shape

    return values, probabilities


import scipy


def calculate_pk_mmk1n(K, N, ro):
    assert(K > 0)
    assert(ro >= 0)
    p = np.zeros((K+1,))

    if (N > 0):
        print(f"\tN packets:{N} buffer size:{K}")
        for k in range(0, K+1):
            p[k] = scipy.special.binom(N, k) * (ro**k) * math.factorial(k)
        p = p / np.sum(p)
        assert(p[0] > 0)  # Normalize to 1
        p = np.clip(p, 0, 1)
        print(f":\tp0:{p[0]:e}")
        print(f"\tpk:{p}")

    return p


def calculate_queue_blocking(K, arrival_rate, mean_service_time,
                             N=None, service_pdf=None, sim_time=None):
    #print("*** Calculate queue blocking ***")
    #print(f"arrival_rate: {arrival_rate:e}")
    #print(f"mean service time: {mean_service_time:e}")
    #assert(mean_service_time > 0)
    #service_rate = 1 / mean_service_time

    ro = arrival_rate * mean_service_time
    print(f"arrival_rate:{arrival_rate:e} mean_service_time:{mean_service_time}")
    """
    if service_pdf is not None:
        # Probability density function was given, calculate parameters for an
        #  M/G/1/K queue
        # TODO: compare results with and without this
        #print(f"Calculating stats for queue of size: {K}")

        #assert(np.sum(service_pdf[1]) == 1.0)
        pi = calculate_pi(K, arrival_rate, service_pdf)

        probability_blocking = 1.0 - (1.0 / (pi[0] + ro))
        #print(f"pi0: {pi[0]} ro:{ro} probability_blocking:{probability_blocking}")

        waiting_sum = 0
        for k in range(1, K):
            waiting_sum += pi[k]
        mean_waiting_time = (waiting_sum/arrival_rate) + ((K/arrival_rate)*(pi[0] + ro - 1)) - mean_service_time
        #print(f"mean waiting time:{mean_waiting_time} mean_service_time:{mean_service_time}")
    else:
    """
    """
    if N is not None:
        print("MM1KN")
        # Calculate p0
        server_utilization = ro  # Rho in queueing theory
        # TODO: we need to reevaluate p[k] for size N-1?
        # TODO: simplify the formula, we don't need to do the factorial.
        #  I think we actually want perm nPr -- or N_Pk
        #  way
        print("Calculating Pk(N-1, 1, K)")
        p = calculate_pk_mmk1n(K, N-1, ro)
        probability_blocking = p[K]
        probability_blocking = np.nan_to_num(probability_blocking)
        #exit()
        print("Calculating Pk(N, 1, K)")
        p = calculate_pk_mmk1n(K, N, ro)
        effective_throughput = arrival_rate * (1 - probability_blocking)  # lambda_e
        #total_length = 0  # L
        queue_length = 0
        for k in range(1, K+1):
            queue_length += (k-1) * p[k]
            #total_length += k * p[k]

        total_length = queue_length
        if total_length > 0:
            total_wait = total_length / effective_throughput  # W
        else:
            total_wait = 0  # W
        mean_waiting_time = total_wait - mean_service_time  # Wq
        effective_utilization = effective_throughput * mean_service_time # rho_e
        #queue_length = total_length - effective_utilization

        print(f"mean service time:{mean_service_time} arrival rate:{arrival_rate:e} server_utilization:{server_utilization}")
        print(f"prob blocking:{probability_blocking:e} total length:{total_length} total waiting time:{total_wait} effective_throughput:{effective_throughput:e}")
        print(f"queue length:{queue_length} estimated waiting time:{mean_waiting_time}")
    """
    #else:
    print("\tMM1K")
    # No probability density function given, assume M/M/K/1 queue
    server_utilization = ro  # Rho in queueing theory
    #print(f"server utilization: {server_utilization}")
    if (server_utilization == 1):
        probability_blocking = 1 / (K+1)
        total_length = K/2
    else:
        probability_blocking = ((1 - ro)*(ro**K)) / (1 - (ro**(K+1)))
        total_length = (server_utilization/(1.0-server_utilization))
        total_length -= (((K+1) * (server_utilization**(K+1))) /
                            (1 - server_utilization**(K+1)))  # L

    #p = np.zeros((K+1))
    #p[0] = (1-ro) / (1-(ro**(K+1)))
    #for k in range(0, K+1):
    #    p[k] = ro**k * p[0]
    #print(f"p:{p}")

    if N is not None and mean_service_time != 0:
        mean_service_rate = (1/mean_service_time)
        print(f"\tservice rate:{mean_service_rate:e}")
        #time_until_saturation = K / (arrival_rate - mean_service_rate)
        #print(f"time until saturation:{time_until_saturation}")
        #time_until_all_sent = N / arrival_rate
        #print(f"N:{N} K:{K} time until all sent:{time_until_all_sent}")

        #if sim_time is not None:
        #steady_state_ratio = (time_until_all_sent - time_until_saturation) / time_until_all_sent
        #steady_state_ratio = np.clip(steady_state_ratio, 0, 1)

        # TODO: this ratio should change depending on the average depth of the link
        #  i.e., the average capacity of downstream links. The number of packets
        #  needs to be greater than this total capacity before we will see any
        #  blocking effects
        # N: Messages sent in flow
        # K: Buffer size
        steady_state_ratio2 = (N - K) / N
        steady_state_ratio2 = np.clip(steady_state_ratio2, 0, 1)
        # TODO: HACK - to disable
        #steady_state_ratio2 = 1.0
        print(f"\tsim time:{sim_time} steady state ratio2:{steady_state_ratio2}")
        #input()
        probability_blocking *= steady_state_ratio2


    # lambda_e == effective throughput
    effective_throughput = arrival_rate * (1 - probability_blocking)
    if total_length > 0:
        total_wait = total_length / effective_throughput
    else:
        total_wait = 0
    mean_waiting_time = total_wait - mean_service_time

    effective_utilization = effective_throughput * mean_service_time # rho_e
    queue_length = total_length - effective_utilization

    probability_blocking = np.clip(probability_blocking, 0.0, 1.0)
    mean_waiting_time = np.clip(mean_waiting_time, 0.0, None)
    queue_length = np.clip(queue_length, 0, None)

    print(f"\tserver_utilization:{server_utilization} total length:{total_length} total wait:{total_wait}")
    print(f"\tprob blocking {probability_blocking} queue length:{queue_length} effective_throughput:{effective_throughput:e} waiting_time:{mean_waiting_time}")
    ##print(f"mean service rate:{mean_service_rate}")

    return probability_blocking, mean_waiting_time, queue_length


def update_buffer_queue(dependencies, link, link_arrival_rates,
                        link_server_time, contention_waiting_time,
                        router_link_flows, path_arrival_rates, flows,
                        prob_link_blocking, link_arrival_count,
                        messages_buffered, max_neuron_processing):
    print("** Update buffer queue **")
    links_out = dependencies.out_edges(link)
    print(f"links_out:{links_out}")
    link_idx = reverse_graph_index(link)

    # Calculate the average server rate of all downstream links
    for downstream in links_out:
        downstream_link = downstream[1]
        print(f"downstream link:{downstream_link}")
        downstream_idx = reverse_graph_index(downstream_link)

        arrival_rate_between_links = 0
        # Find all flows through this link, going to the downstream link
        flows_in_both_links = list(set(router_link_flows[link]) &
                               set(router_link_flows[downstream_link]))
        for f in flows_in_both_links:
            path_idx = flows[f]
            arrival_rate_between_links += path_arrival_rates[path_idx[0],
                                                             path_idx[1]]
            print(f"flow:{f} path idx:{path_idx} arrival:{path_arrival_rates[path_idx[0], path_idx[1]]:e}")

        print(f"contention waiting time:{contention_waiting_time[downstream_idx]:e}")
        print(f"arrival rate between links {link_idx}->{downstream_idx} = "
              f"{arrival_rate_between_links:e}")
        print(f"total link arrival rate: {link_arrival_rates[link_idx]:e}")
        print(f"scaled contention time: {((contention_waiting_time[downstream_idx]*arrival_rate_between_links) / link_arrival_rates[link_idx]):e}")
        print(f"prob blocking:{prob_link_blocking[downstream_idx]}")
        #add = (contention_waiting_time[downstream_idx] + \
        #    (1 / (1 - prob_link_blocking[downstream_idx]))) * \
        #    arrival_rate_between_links
        #link_server_time[link_idx] += \
        #    (contention_waiting_time[downstream_idx] + \
        #    (1 / (1 - prob_link_blocking[downstream_idx]))) * \
        #    arrival_rate_between_links
        link_server_time[link_idx] += \
            (arrival_rate_between_links) * contention_waiting_time[downstream_idx]

    # Normalize server time with respect with the total flow going through this
    # link
    if (len(links_out) > 0) and (link_arrival_rates[link_idx] > 0):
        link_server_time[link_idx] /= link_arrival_rates[link_idx]
    print(f"buffer server time: {link_server_time[link_idx]}")
    assert(link_server_time[link_idx] != math.nan)
    assert(link_server_time[link_idx] < 1.0)
    # else this is a leaf link, and the service time is already given by the path

    #"""
    if link_idx[2] == 0 or link_idx[2] == 2:  # Network tile to network tile
        buffer_size = 16
    elif link_idx[2] == 1 or link_idx[2] == 3:
        buffer_size = 10
    elif (link_idx[2] >= 4 and link_idx[2] < 8):  # Cores to network
        buffer_size = 8
    elif link_idx[2] >= 8:  # Network to cores
        buffer_size = 24
    #"""

    prob_blocked, link_waiting_time, queue_length = \
        calculate_queue_blocking(buffer_size, link_arrival_rates[link_idx],
                                 link_server_time[link_idx],
                                 N=link_arrival_count[link_idx],
                                 service_pdf=None,
                                 sim_time=max_neuron_processing)
    #prob_blocked, link_waiting_time, queue_length = \
    #    calculate_queue_blocking(buffer_size, link_arrival_rates[link_idx],
    #                             link_server_time[link_idx],
    #                             N=None,
    #                             service_pdf=None,
    #                             sim_time=max_neuron_processing)

    """
    prob_blocked1, link_waiting_time1 = \
        calculate_queue_blocking(10, link_arrival_rates[link_idx],
                                 link_server_time[link_idx],
                                 N=link_arrival_count[link_idx],
                                 service_pdf=None)
    prob_blocked2, link_waiting_time2 = \
        calculate_queue_blocking(100, link_arrival_rates[link_idx],
                                 link_server_time[link_idx],
                                 N=link_arrival_count[link_idx],
                                 service_pdf=None)
    print(f"prob block 1:{prob_blocked1} 2:{prob_blocked2}")
    """
    #input()

    print(f"link buffer waiting:{link_waiting_time:e}")
    prob_link_blocking[link_idx] = np.clip(prob_blocked, 0, 1)
    messages_buffered[link_idx] = np.clip(queue_length, 0, None)
    return link_waiting_time


def update_contention_queue(dependencies, link, link_arrival_rates,
                            contention_waiting_time,
                            prob_link_blocking,
                            link_server_time,
                            sim_time):
    print("** Update contention queue **")
    links_in = dependencies.in_edges(link)
    print(f"links_in:{links_in}")
    link_in_count = len(links_in)
    print(f"link in count:{link_in_count}")
    link_idx = reverse_graph_index(link)

    # The server time of the queue is just the average time that the link
    #  is blocked for
    time_blocked = 0
    if (link_in_count >= 1) and link_arrival_rates[link_idx] > 0:
        # This equation doesn't work / make sense
        #time_blocked = (1 / link_arrival_rates[link_idx]) * (1 / (1 - prob_link_blocking[link_idx]))
        time_blocked = prob_link_blocking[link_idx] * link_server_time[link_idx]
        print(f"time_blocked:{time_blocked} arrival:{link_arrival_rates[link_idx]:e} prob_block:{prob_link_blocking[link_idx]}")
        assert(time_blocked != math.nan)
        assert(time_blocked >= 0)
        contention_server_time = time_blocked
        print(f"prob blocked: {prob_link_blocking[link_idx]} time blocked:{contention_server_time:e}")

        _, contention_waiting_time[link_idx], _ = \
            calculate_queue_blocking(link_in_count, link_arrival_rates[link_idx],
                                    contention_server_time, None,
                                    sim_time=sim_time)
    print(f"contention waiting time:{contention_waiting_time[link_idx]}")
    print ("**end of contention**")

    return


n_tiles = 32
router_link_names = ("north", "east", "south", "west",
                    "core_1_to_net", "core_2_to_net",
                    "core_3_to_net", "core_4_to_net",
                    "net_to_core_1", "net_to_core_2",
                    "net_to_core_3", "net_to_core_4")


def sim_delay_hops(flows):
    """Calculate delay of each message based on the hop count

    The simplest delay calculation scheme, beyond just setting network delay=0
    """
    delays = np.zeros((128, 128))
    cost_per_hop = 4.1e-9  # s
    for src_core, dest_core in flows:
        src_tile = src_core // 4
        dest_tile = dest_core // 4
        src_x = src_tile // 4
        src_y = src_tile % 4
        dest_x = dest_tile // 4
        dest_y = dest_tile % 4
        x_hops = abs(src_x - dest_x)
        y_hops = abs(src_y - dest_y)
        delays[src_core][dest_core] = (x_hops + y_hops) * cost_per_hop

    return delays


def sim_schedule_event_based(messages):
    # Model the queues in each router link as a set of FIFOs. Here we will
    #  essentially schedule every single event as messages go from router to
    #  router.
    noc_buffers = []
    last_sent_buffers = np.zeros((8, 4, 12))
    pending_messages = [None for _ in range(0, 128)]
    receiving_messages = [None for _ in range(0, 128)]
    processed_messages = [[] for _ in range(0, 128)]
    all_messages = set()

    total_messages_sent = 0
    for m_q in messages:
        m_q.reverse()
        total_messages_sent += len(m_q)
        for m in m_q:
            all_messages.add(m)

    for _ in range(0, 8):
        y_buffers = []
        for _ in range(0, 4):
            link_buffers = [[] for _ in range(0, 12)]
            y_buffers.append(link_buffers)
        noc_buffers.append(y_buffers)

    class State(enum.Enum):
        SEND = 1
        RECEIVE = 2

    def update_noc(t):
        # TODO: this is probably simpler for handling pending messages
        #  Go through all cores and check for any pending messages
        # Go through all pending messages and router links and check if any
        #  messages can now progress further through the NoC
        for x in range(0, 8):
            for y in range(0, 4):
                for link in range(4, 8):
                    if len(noc_buffers[x][y][link]):
                        src_core = (x*16)+(y*4)+(link%4)
                        m = pending_messages[src_core]
                        if m is not None:
                            if len(noc_buffers[x][y][link]) < 8:
                                m.pending = False
                                pending_messages[src_core] = None
                                send_message(m, t)
                                # Schedule the next messages if applicable
                                if len(messages[m.src_core]) > 0:
                                    next_message = messages[m.src_core].pop()
                                    logging.debug(f"After pending message handled, scheduled m:{next_message}")
                                    message_sent = send_message(next_message, t)
                                    logging.info(f"Message m:{next_message} sent:{message_sent} at t={t}")
                                    if next_message.pending:
                                        pending_messages[next_message.src_core] = next_message
                                    elif len(messages[next_message.src_core]) > 0:
                                        logging.debug("Message sent ok, schedule next message in queue")
                                        # TODO: when we schedule the next message,
                                        #  that will occur some time after the prev
                                        #  message.
                                        schedule_next_message(next_message, t)
                #for link in (0, 1, 2, 3, 8, 9, 10, 11):
                #    if len(noc_buffers[x][y][link]) > 0:
                #        update_message(noc_buffers[x][y][link][-1], t)

    def sanity_check_noc():
        logging.trace("Sanity checking NoC state")
        for src in range(0, 128):
            # Assert that we can't have a pending message and space in that
            #  buffer
            tile = src // 4
            x = tile // 4
            y = tile % 4
            link = 4 + (src % 4)
            if pending_messages[src]:
                buffer_size = 8
                logging.trace(f"core:{src} link:{link} "
                             f"pending:{pending_messages[src]} "
                             f"buf out len:{len(noc_buffers[x][y][link])}")
                assert(len(noc_buffers[x][y][link]) >= buffer_size)

        logging.trace("Checking to see links aren't jammed")
        for x in range(0, 8):
            for y in range(0, 4):
                for link in range(0, 4):
                    if len(noc_buffers[x][y][link]) > 0:
                        m = noc_buffers[x][y][link][-1]
                        logging.info(f"checking m:{m}")
                        dest_tile = m.dest_core // 4
                        dest_x = dest_tile // 4
                        dest_y = dest_tile % 4

                        if dest_y > y:
                            next_x, next_y = x, y+1
                        elif dest_x > x:
                            next_x, next_y = x+1, y
                        elif dest_y < y:
                            next_x, next_y = x, y-1
                        elif dest_x < x:
                            next_x, next_y = x-1, y
                        else:
                            next_x, next_y = x, y

                        if x < dest_x:
                            next_link = 1  # East
                        elif x > dest_x:
                            next_link = 3  # West
                        elif y < dest_y:
                            next_link = 0  # North
                        elif y > dest_y:
                            next_link = 2  # South
                        else:
                            next_link = 8 + (m.dest_core % 4)

                        buffer_sizes = (16, 10, 16, 10, 8, 8, 8, 8, 24, 24, 24, 24)
                        logging.trace(f"link {x,y,link}->{next_x,next_y,next_link} "
                                     f"buffer len:{len(noc_buffers[next_x][next_y][next_link])} "
                                     f"max len:{buffer_sizes[next_link]}")
                        assert(len(noc_buffers[next_x][next_y][next_link]) == \
                               buffer_sizes[next_link])
                for link in range(8, 12):
                    for m in noc_buffers[x][y][link]:
                        if (m.dest_core % 4) != (link % 4):
                            logging.trace(f"m:{m} and link:{link}")
                        assert((m.dest_core % 4) == (link % 4))
        return

    def get_next_downstream_link(m):
        x, y, link = m.pos
        if link is None:
            # if pos is None, the next link is the first one
            src_tile = m.src_core // 4
            next_x = src_tile // 4
            next_y = src_tile % 4
            next_link = 4 + (m.src_core % 4)
            logging.debug(f"m:{m} entered at link {next_x, next_y, next_link}")
        else:
            assert(link < 8)  # Already at destination link
            dest_tile = m.dest_core // 4
            dest_x = dest_tile // 4
            dest_y = dest_tile % 4
            if x < dest_x:
                next_x, next_y = x+1, y
            elif x > dest_x:
                next_x, next_y = x-1, y
            elif y < dest_y:
                next_x, next_y = x, y+1
            elif y > dest_y:
                next_x, next_y = x, y-1
            else:
                next_x, next_y = x, y

            # Work out the direction of the next link, assuming XY routing
            if x < dest_x:
                next_link = 1  # East
            elif x > dest_x:
                next_link = 3  # West
            elif y < dest_y:
                next_link = 0  # North
            elif y > dest_y:
                next_link = 2  # South
            else:
                next_link = 8 + (m.dest_core % 4)

        # If the next link is full, it will block any sending link
        buffer_sizes = (16, 10, 16, 10, 8, 8, 8, 8, 24, 24, 24, 24)
        if len(noc_buffers[next_x][next_y][next_link]) >= buffer_sizes[next_link]:
            logging.debug(f"next link {next_x,next_y,next_link} is full i.e., blocking")
            return None, None, None  # To indicate the next link is full
        else:
            logging.debug(f"next link: {next_x, next_y, next_link}")
            return next_x, next_y, next_link


    def receive_message(m, t):
        # Schedule an event to update the NoC once we finish processing
        logging.debug(f"creating a RECEIVE event for m:{m} at t={t}")
        E = {"message": m, "type": State.RECEIVE}
        priority.append((E, t + m.receive_delay))
        priority.sort(key=lambda x: x[1], reverse=True)
        assert(receiving_messages[m.dest_core] is None)
        receiving_messages[m.dest_core] = m
        logging.trace(f"receiving messages:{receiving_messages}")
        return

    def send_message(m, t):
        logging.debug(f"Sending message {m}@t:{t}")
        dest_pos = get_next_downstream_link(m)
        dest_x, dest_y, dest_link = dest_pos

        if dest_link is None:
            logging.debug(f"First link is busy, set message m:{m} as pending")
            # Message is being blocked by NoC, we need to push this
            m.pending = True
            pending_messages[m.src_core] = m
            logging.debug(f"Pending:{pending_messages}")
            return False
        else:
            m.pending = False

        m.pos = dest_pos
        logging.debug(f"m:{m} m.pos:{m.pos}")
        x, y, link = m.pos
        while (dest_link is not None and
               (len(noc_buffers[dest_x][dest_y][dest_link]) == 0) and
               dest_link < 8):
            # Update the router tracking of which links have been used
            logging.debug(f"Router link ({x},{y},{link}) "
                         f"accessed router @ t={t}")
            last_sent_buffers[x, y, link] = t
            dest_pos = get_next_downstream_link(m)
            dest_x, dest_y, dest_link = dest_pos

            if dest_link is not None:
                m.pos = dest_pos
                logging.debug(f"m:{m} m.pos:{m.pos}")
                x, y, link = m.pos

        if (len(noc_buffers[x][y][link]) == 0
            and receiving_messages[m.dest_core] is None and (link >= 8)):
            # We successfully sent the message across the NoC, no blocking
            #  We can start receiving the message
            receive_message(m, t)
            logging.info(f"Message m:{m} sent and received at core {m.dest_core}")
        else:
            logging.info(f"Message m:{m} sent")
            logging.debug(f"len buffer[{x,y,link}]:{len(noc_buffers[x][y][link])}")
            noc_buffers[x][y][link].insert(0, m)

        logging.debug(f"Message m:{m} stored at pos q{m.pos}")
        return True

    def update_message(m, t):
        # Update a message in the NoC. Assuming that a queue upstream has become
        #  available, this means the message can possibly make its way to its
        #  destination core now
        src_x, src_y, src_link = m.pos
        hops = 0

        if src_link < 8:
            logging.debug(f"Get next link for m:{m}")
            # If we are not at a dest link, so somewhere inside the NoC, find
            #  the next link for message m
            dest_pos = get_next_downstream_link(m)
            dest_x, dest_y, dest_link = dest_pos

            logging.debug(f"Updating m:{m} at link {m.pos}")
            if dest_link is None:
                logging.debug("Next link is blocked, message can't go anywhere")
                return False
        else:
            # Message m is at its destination link waiting for the dest core
            #  to become available
            logging.debug(f"m:{m} already at dest link")
            dest_x, dest_y, dest_link = src_x, src_y, src_link

        assert(dest_x is not None)
        assert(dest_y is not None)
        assert(dest_link is not None)

        # While there are free spaces downstream, this message will hop across
        #  the NoC
        while (dest_link is not None and
               len(noc_buffers[dest_x][dest_y][dest_link]) == 0 and
               dest_link < 8):
            logging.debug(f"m:{m} found buffer space at link {dest_x,dest_y,dest_link}")
            m.pos = dest_pos
            dest_pos = get_next_downstream_link(m)
            assert(dest_pos != m.pos)
            dest_x, dest_y, dest_link = dest_pos
            x, y, link = m.pos
            last_sent_buffers[x, y, link] = t
            hops += 1

        if dest_link is not None:
            m.pos = dest_x, dest_y, dest_link

        if dest_link is not None and dest_link >= 8:
            # Message made it to its dest link
            logging.debug(f"dest:{dest_link} "
                         f"len:{len(noc_buffers[dest_x][dest_y][dest_link])} "
                         f"rx:{receiving_messages[m.dest_core]}")

            if (len(noc_buffers[dest_x][dest_y][dest_link]) == 0 and
                receiving_messages[m.dest_core] is None):
                # We successfully sent the message across the NoC, and the
                #  receiving core is free to receive the new message
                logging.info(f"message m:{m} started being received at core:{m.dest_core} at t={t}")
                receive_message(m, t)
                # Don't store the message in any buffer or queue, it was
                #  received and processed
                m.pos = None, None, None
                dest_x, dest_y, dest_link = None, None, None
            else:
                logging.debug(f"message waiting to be received, don't receive it")

        noc_buffers[src_x][src_y][src_link].pop()
        x, y, link = m.pos
        if link is not None:
            noc_buffers[x][y][link].insert(0, m)

        logging.debug(f"message m:{m} moved from {(src_x,src_y,src_link)}->"
                        f"{(x,y,link)} dest:{dest_x,dest_y,dest_link}")
        logging.debug(f"message m:{m} moved {hops} hops")

        return True

    def schedule_next_message(m, t):
        # Schedule the next message
        assert(len(messages[m.src_core]) > 0)
        next_message = messages[m.src_core][-1]
        assert(next_message != m)

        logging.debug(f"Scheduling SEND message for m:{next_message}, from core:{m.src_core}")
        E = {"message": next_message, "type": State.SEND}
        priority.append((E, t + next_message.generation_delay))
        priority.sort(key=lambda x: x[1], reverse=True)
        return

    def get_next_upstream_link(x, y, link):
        """Figure out the upstream link to grab a message from"""

        # First figure out the src router position
        if link < 8 and link >= 4:  # This is a link from src core->NoC
            # There are no more upstream links
            return None, None, None
        elif link >= 8:  # This is a link from NoC->dest core
            # The same router arbitrates the connecting link
            next_x, next_y = x, y
        elif link == 0:  # North
            next_x, next_y = x, y-1
        elif link == 1:  # East
            next_x, next_y = x-1, y
        elif link == 2:  # South
            next_x, next_y = x, y+1
        elif link == 3:  # West
            next_x, next_y = x+1, y

        logging.debug(f"last sent:{last_sent_buffers[next_x,next_y,:]}")
        masked = np.ones((8,)) * np.inf
        arbitration_count = 0
        # Go to every link and figure out if the next message is contending over
        #  the current link
        for buf in noc_buffers[next_x][next_y][:8]:
            if len(buf) > 0:
                m = buf[-1]
                logging.debug(f"checking contention m:{m} pos:{x,y,link}")
                target_x, target_y, target_link = get_next_downstream_link(m)
                logging.debug(f"m:{m} targeting link {target_x,target_y,target_link}")
                if (target_x == x) and (target_y == y) and (target_link == link):
                    logging.debug(f"m:{m} competing over the dest link, buffer len:{len(noc_buffers[target_x][target_y][target_link])}")
                    masked[m.pos[2]] = \
                        last_sent_buffers[m.pos[0], m.pos[1], m.pos[2]]
                    arbitration_count += 1

        logging.debug(f"links competing:{arbitration_count} over buffer len:{len(noc_buffers[x][y][link])}")
        logging.debug(f"masked last sent times:{masked}")

        if arbitration_count > 0:
            next_link = np.argmin(masked)
            logging.debug(f"Chosen link:{np.argmin(masked)} with t={masked[next_link]}")
            logging.debug(f"For link[{x}][{y}][{link}] "
                        f"next arbitrated link is [{next_x}][{next_y}][{next_link}]")
            return next_x, next_y, next_link
        else:
            logging.debug("No competing links upstream")
            return None, None, None

    def unblock_link(x, y, link, t):
        # The link was blocked until we updated it
        logging.debug(f"Link {x, y, link} was blocking and is now unblocked")
        # Find the next link and get its first message, if it exists
        if link >= 4 and link < 8:
            src_core = (x*16)+(y*4)+(link%4)
            logging.debug(f"src core {src_core} unblocked, checking pending")
            m = pending_messages[src_core]
            pending_messages[src_core] = None
            logging.debug(f"Pending message:{m}")
            next_x, next_y, next_link = None, None, None

            if m is not None:  # Pending message found that will now be sent
                logging.info(f"Pending message m:{m} entered NoC ")
                m.pending = False
                # Message has been waiting as pending, so send it out
                #  immediately
                send_message(m, t)
                # Schedule the next messages if applicable, or put it into the
                #  pending state
                if len(messages[m.src_core]) > 0:
                    next_message = messages[m.src_core].pop()
                    logging.debug(f"After pending message handled, scheduled m:{next_message}")
                    message_sent = send_message(next_message, t)
                    logging.info(f"Message m:{next_message} sent:{message_sent} at t={t}")
                    if not next_message.pending and len(messages[next_message.src_core]) > 0:
                        logging.debug("Message sent ok, schedule next message in queue")
                        schedule_next_message(next_message, t)
        else:  # A router->router link
            buffered_messages = buffer_sizes[link]
            while len(noc_buffers[x][y][link]) < buffer_sizes[link]:
                logging.debug(f"Getting upstream link for {x,y,link}")
                next_pos = get_next_upstream_link(x, y, link)
                logging.debug(f"Upstream link is:{next_pos}")
                # Get the next link to win arbitration to send to the current
                #  link
                next_x, next_y, next_link = next_pos

                if (next_link is not None and len(noc_buffers[next_x][next_y][next_link]) > 0):
                    logging.debug(f"Checking upstream link messages for link {next_x,next_y,next_link}")
                    buffered_messages = len(noc_buffers[next_x][next_y][next_link])

                    upstream_link_sending = True
                    messages_sent = 0
                    # Unblock as messages in the buffer queue as possible
                    while (upstream_link_sending and
                        len(noc_buffers[next_x][next_y][next_link]) > 0):
                        assert(noc_buffers[next_x][next_y][next_link][-1].pos == \
                                (next_x, next_y, next_link))
                        upstream_link_sending = \
                            update_message(noc_buffers[next_x][next_y][next_link][-1], t)
                        logging.debug(f"Upstream link ok to send more messages:{upstream_link_sending}")
                        messages_sent += 1

                    logging.debug(f"{messages_sent} messages sent after unblocking link {next_x, next_y, next_link}")
                    logging.debug(f"buffered messages:{buffered_messages} buffer_sizes[]={buffer_sizes[next_link]}")
                    if buffered_messages >= buffer_sizes[next_link]:
                        unblock_link(next_x, next_y, next_link, t)
                    else:
                        logging.debug(f"Upstream link {next_x, next_y, next_link} "
                                    f"{buffered_messages} messages remaining before updating, "
                                    f"{len(noc_buffers[next_x][next_y][next_link])} after updating")
                        break
                    x, y, link = next_x, next_y, next_link
                else:  # There is no upstream link with messages
                    src_core = (x*16)+(y*4)+(link%4)
                    logging.debug(f"src core {src_core} checking pending")
                    m = pending_messages[src_core]
                    if m is not None:
                        logging.info(f"Pending message m:{m} entered NoC ")
                        m.pending = False
                        send_message(m, t)
                        # Schedule the next messages if applicable
                        if len(messages[m.src_core]) > 0 :
                            next_message = messages[m.src_core].pop()
                            logging.debug(f"After pending message handled, scheduled m:{next_message}")
                            message_sent = send_message(next_message, t)
                            logging.info(f"Message m:{next_message} sent:{message_sent} at t={t}")
                            if not next_message.pending and len(messages[next_message.src_core]) > 0:
                                schedule_next_message(next_message, t)
                    break

        return

    # Initialize priority queue by pushing events to send the first messages
    priority = []
    for core in range(0, 128):
        if len(messages[core]) > 0:
            first_message = messages[core][-1]
            E = {"message": first_message, "type": State.SEND}
            priority.append((E, first_message.generation_delay))
    priority.sort(key=lambda x: x[1], reverse=True)
    last_t = 0

    while len(priority) > 0:
        E, t = priority.pop()
        logging.info("")
        logging.info(f"************ EVENT E:{E} @ t={t} ************")
        m = E["message"]

        if E["type"] == State.SEND:
            # First, figure out if we are able to send a message or not. If the
            #  NoC is blocking, put the message into a pending state. TODO:
            #  I thin the receive event will need a pending buffer per core
            #  to know.
            m = messages[m.src_core].pop()
            assert(m == E["message"])
            logging.debug(f"Sending message m:{m} pending:{m.pending}")

            send_message(m, t)

            logging.debug(f"m:{m} pending:{m.pending}")
            logging.trace(f"pending:{pending_messages}")
            # Message has been sent successfully, schedule the next message
            if m.pending:
                pending_messages[m.src_core] = m
            elif len(messages[m.src_core]) > 0:
                logging.debug("Message sent ok, schedule next message in queue")
                schedule_next_message(m, t)

        elif E["type"] == State.RECEIVE:  # State.Receive
            logging.info(f"Message m:{m} finished being received by core:{m.dest_core} at t={t}")
            receiving_messages[m.dest_core] = None
            m.t_finished_processing = t
            processed_messages[m.dest_core].append(m)
            dest_tile = m.dest_core // 4
            x = dest_tile // 4
            y = dest_tile % 4
            link = 8 + (m.dest_core % 4)
            logging.debug(f"Checking link[{x,y,link} for messages, "
                          f" count={len(noc_buffers[x][y][link])}")
            logging.trace(f"link: {noc_buffers[x][y][link]}")
            if len(noc_buffers[x][y][link]) > 0:  # Another message in the queue
                logging.debug("link is not empty, getting next message")
                # Start processing the next message in the receive queue
                buffered_messages = len(noc_buffers[x][y][link])
                next_m = noc_buffers[x][y][link].pop()
                logging.debug(f"next m:{next_m}")
                assert(next_m.pos == (x, y, link))
                logging.info(f"message m:{next_m} started being received at core:{m.dest_core} at t={t}")
                receive_message(next_m, t)

                # While router links are unblocked, pop one element and look
                #  downstream for the next message(s) to fill the gap
                buffer_sizes = (16, 10, 16, 10, 8, 8, 8, 8, 24, 24, 24, 24)
                if buffered_messages >= buffer_sizes[link]:
                    unblock_link(x, y, link, t)
            else:
                logging.debug("Next link is empty, finished receiving")
        else:
            raise Exception(f"Error: Event:{E} not supported")
        update_noc(t)
        assert(t >= last_t)
        last_t = max(t, last_t)
        # TODO: sanity check is broken or overly sensitive in some casses

        # Check here that the total messages in the system are consistent
        total_in_system = 0
        check_messages = set()
        for m_q in messages:
            total_in_system += len(m_q)
            for m in m_q:
                check_messages.add(m)
        for x in range(0, 8):
            for y in range(0, 4):
                for link in range(0, 12):
                    total_in_system += len(noc_buffers[x][y][link])
                    for m in noc_buffers[x][y][link]:
                        check_messages.add(m)
        for m_q in processed_messages:
            total_in_system += len(m_q)
            for m in m_q:
                check_messages.add(m)
            check_messages.add(m)
        for m in receiving_messages:
            if m is not None:
                total_in_system += 1
                check_messages.add(m)
        for m in pending_messages:
            if m is not None:
                total_in_system += 1
                check_messages.add(m)

        logging.debug(f"total in system:{total_in_system} total messages:{total_messages_sent}")
        missing = set(all_messages) ^ set(check_messages)
        logging.debug(f"missing messages: {missing}")
        #assert(total_in_system == total_messages_sent)
        #sanity_check_noc()
        # ** end of event loop **

    # Do final sanity checks
    for c, m in enumerate(messages):
        logging.trace(f"{c}:{m}")

    for x in range(0, 8):
        for y in range(0, 4):
            for link in range(0, 4):
                logging.trace(f"{x,y,link}: {noc_buffers[x][y][link]}")

    logging.trace(f"receiving:{receiving_messages}")
    logging.trace(f"pending:{pending_messages}")
    total_messages_received = 0
    """
    for c, m_q in enumerate(processed_messages):
        print(f"Dest core {c}:")
        for m in m_q:
            total_messages_received += 1
            print(f"\t\tm:{m} gen delay:{m.generation_delay} rcv delay:{m.receive_delay} processed:{m.t_finished_processing}")
    """
    print(f"total received:{total_messages_received} total sent:{total_messages_sent}")

    #assert(all([p is None for p in pending_messages]))
    #assert(all([r is None for r in receiving_messages]))
    #assert(all([len(m) == 0 for m in messages]))
    #assert(total_messages_received == total_messages_sent)

    print("*** Finished modeling NoC ***")
    #exit()
    return last_t


def sim_delay_mm1k(df):
    n_cores = 128
    path_counts = np.zeros((n_cores, n_cores), dtype=int)  # [src core, dest core]
    path_arrival_latencies = np.zeros((n_cores, n_cores))
    path_server_mean_latencies = np.zeros((n_cores, n_cores))
    router_link_flows = [[] for _ in range(0, 8*4*12)]

    path_server_latencies = [[] for _ in range(0, 128)]
    path_server_pdf = [[] for _ in range(0, 128)]
    for s in range(0, 128):
        path_server_latencies[s] = [[] for _ in range(0, 128)]
        path_server_pdf[s] = [[] for _ in range(0, 128)]

    router_link_counts = np.zeros((8, 4, 12))
    router_link_arrival_rates = np.zeros((8, 4, 12))

    dependencies = nx.DiGraph()
    # Flatten graph index into 1-dimension
    dependencies.add_nodes_from(range(0, 8*4*12))

    # Parse the rate of messages between two cores in the network and build a
    #  dependency graph
    for _, row in df.iterrows():
        src_tile, src_core = hw_str_to_core(row["src_hw"])
        dest_tile, dest_core = hw_str_to_core(row["dest_hw"])
        path_arrival_latencies[src_core, dest_core] += row["generation_delay"]
        #if row["processing_latency" != 0]:
        #    print(f"latency:{path_arrival_latencies[src_core, dest_core]}")
        assert(row["processing_latency"] >= 0)
        assert(row["processing_latency"] != math.nan)
        path_server_mean_latencies[src_core, dest_core] += \
            row["processing_latency"]
        assert(np.all(path_server_mean_latencies >= 0))
        path_counts[src_core, dest_core] += 1
        path_server_latencies[src_core][dest_core].append(row["processing_latency"])

    #assert(np.all(path_server_mean_latencies >= 0))
    #assert(not np.any(np.isnan(path_server_mean_latencies)))
    #assert(not np.any(np.isnan(path_counts)))
    #assert(np.all(path_counts >= 0))

    path_server_mean_latencies = np.divide(path_server_mean_latencies,
                                           path_counts, where=(path_counts>=1),
                                           dtype=float)
    np.clip(path_server_mean_latencies, 0.0, None)
    path_server_mean_latencies = np.nan_to_num(path_server_mean_latencies)

    #path_arrival_rates = np.divide(path_counts, np.max(path_arrival_latencies),
    #                               where=path_counts>0)  # packets/s
    path_arrival_rates = np.divide(path_counts, path_arrival_latencies*10,
                                   where=path_counts>0)  # packets/s

    # Treat the neuron processing time as the simulation window
    max_neuron_processing = np.max(np.sum(path_arrival_latencies, axis=1))

    #print(np.ma.masked_less(path_server_rates, 1.0))
    flows = np.argwhere(path_counts >= 1)
    # Figure out the distribution of latencies for each flow
    for s in range(0, 128):
        for d in range(0, 128):
            path_server_pdf[s][d] = create_pdf(path_server_latencies[s][d])

    #print(path_rates)

    # *** Router link model variables ***
    prob_link_blocking = np.zeros((8, 4, 12))
    mean_link_waiting_time = np.zeros((8, 4, 12))
    mean_link_service_time = np.zeros((8, 4, 12))
    messages_buffered = np.zeros((8, 4, 12))

    # *** Contention model variables ***
    contention_waiting_time = np.zeros((8, 4, 12))

    # Now parse all the links between two cores. Every router has 12 links in
    #  the NoC: 4 links in the going North, East, south and West. Every router
    #  has 4 links going from the core to the router, and 4 links going out to
    #  the cores. Compute the packet arrival rates for all links on all paths
    for idx, flow in enumerate(flows):
        src_core, dest_core = flow
        assert(src_core < 128)
        assert(dest_core < 128)

        src_tile = src_core // 4
        dest_tile = dest_core // 4
        assert(src_tile < 32)
        assert(dest_tile < 32)

        src_x = src_tile // 4
        src_y = src_tile % 4
        dest_x = dest_tile // 4
        dest_y = dest_tile % 4
        path_rate = path_arrival_rates[src_core, dest_core]
        path_count = path_counts[src_core, dest_core]

        # Account for all the hops in between the src and dest tile
        initial_direction = 4 + (src_core % 4)
        router_link_counts[src_x, src_y, initial_direction] += path_count
        router_link_arrival_rates[src_x, src_y, initial_direction] += path_rate
        router_link_flows[graph_index(src_x, src_y, initial_direction)].append(idx)
        prev_link = (src_x, src_y, initial_direction)

        while src_x < dest_x:
            src_x += 1
            router_link_counts[src_x, src_y, 1] += path_count  # east
            router_link_arrival_rates[src_x, src_y, 1] += path_rate
            router_link_flows[graph_index(src_x, src_y, 1)].append(idx)
            dependencies.add_edge(graph_index(*prev_link),
                                graph_index(src_x, src_y, 1))
            prev_link = (src_x, src_y, 1)

        while src_x > dest_x:
            src_x -= 1
            router_link_counts[src_x, src_y, 3] += path_count  # west
            router_link_arrival_rates[src_x, src_y, 3] += path_rate
            router_link_flows[graph_index(src_x, src_y, 3)].append(idx)
            dependencies.add_edge(graph_index(*prev_link),
                                graph_index(src_x, src_y, 3))
            prev_link = (src_x, src_y, 3)

        while src_y < dest_y:
            src_y += 1
            router_link_counts[src_x, src_y, 0] += path_count  # north
            router_link_arrival_rates[src_x, src_y, 0] += path_rate
            router_link_flows[graph_index(src_x, src_y, 0)].append(idx)
            dependencies.add_edge(graph_index(*prev_link),
                                graph_index(src_x, src_y, 0))
            prev_link = (src_x, src_y, 0)

        while src_y > dest_y:
            src_y -= 1
            router_link_counts[src_x, src_y, 2] += path_count  # south
            router_link_arrival_rates[src_x, src_y, 2] += path_rate
            router_link_flows[graph_index(src_x, src_y, 2)].append(idx)
            dependencies.add_edge(graph_index(*prev_link),
                                graph_index(src_x, src_y, 2))
            prev_link = (src_x, src_y, 2)

        final_direction = 8 + (dest_core % 4)
        router_link_counts[src_x, src_y, final_direction] += path_count
        router_link_arrival_rates[src_x, src_y, final_direction] += path_rate
        router_link_flows[graph_index(src_x, src_y, final_direction)].append(idx)
        dependencies.add_edge(graph_index(*prev_link),
                            graph_index(src_x, src_y, final_direction))

        # Track the service time at the receiving link
        mean_link_service_time[src_x, src_y, final_direction] = \
            path_server_mean_latencies[src_core, dest_core]
        #print(mean_contention_delay[src_x, src_y, final_direction])


    # Create dependency graph of all links
    sorted_links = list(nx.topological_sort(dependencies))
    sorted_link_info = [reverse_graph_index(n) for n in sorted_links]
    nx.nx_agraph.write_dot(dependencies, "runs/dependencies.dot")
    D = nx.nx_agraph.to_agraph(dependencies)
    # If we want to plot dependencies
    #plt.figure()
    #pos=graphviz_layout(dependencies, prog='dot')
    #nx.draw(dependencies, pos, with_labels=False, arrows=True, node_size=0.5)


    # Now iterate over all flows again in reverse topological order, and
    #  calculate the contention delay followed by the service time at each
    #  link queue. Calculate the total arrival rate at each queue.
        # Calculate the average contention delay for this link


    sorted_links.reverse()
    #print(sorted_links)
    for link in sorted_links:
        print(f"UPDATING LINK: {link}")
        link_wait_time = \
            update_buffer_queue(dependencies, link, router_link_arrival_rates,
                                mean_link_service_time,
                                contention_waiting_time,
                                router_link_flows,
                                path_arrival_rates, flows,
                                prob_link_blocking,
                                router_link_counts,
                                messages_buffered,
                                max_neuron_processing)
        mean_link_waiting_time[reverse_graph_index(link)] = link_wait_time
        update_contention_queue(dependencies, link, router_link_arrival_rates,
                                contention_waiting_time, prob_link_blocking,
                                mean_link_service_time,
                                max_neuron_processing)


    # Now go through all flows, and calculate the delay for each flows path
    flow_latencies = np.zeros((128, 128))
    for flow in flows:
        src_core, dest_core = flow
        src_tile = src_core // 4
        dest_tile = dest_core // 4
        assert(src_tile < 32)
        assert(dest_tile < 32)

        src_x = src_tile // 4
        src_y = src_tile % 4
        dest_x = dest_tile // 4
        dest_y = dest_tile % 4

            # Account for all the hops in between the src and dest tile
        initial_direction = 4 + (src_core % 4)
        flow_latencies[src_core, dest_core] += \
            mean_link_waiting_time[src_x, src_y, initial_direction]

        while src_x < dest_x:
            src_x += 1
            flow_latencies[src_core, dest_core] += \
                mean_link_waiting_time[src_x, src_y, 1]

        while src_x > dest_x:
            src_x -= 1
            flow_latencies[src_core, dest_core] += \
                mean_link_waiting_time[src_x, src_y, 3]

        while src_y < dest_y:
            src_y += 1
            flow_latencies[src_core, dest_core] += \
                mean_link_waiting_time[src_x, src_y, 0]

        while src_y > dest_y:
            src_y -= 1
            flow_latencies[src_core, dest_core] += \
                mean_link_waiting_time[src_x, src_y, 2]

        final_direction = 8 + (dest_core % 4)
        flow_latencies[src_core, dest_core] += \
            mean_link_waiting_time[dest_x, dest_y, final_direction]

    """
    # Plot a heat map
    plt.figure()
    plt.title("Path Arrival Rates")
    plt.xlabel("Source Core")
    plt.ylabel("Destination Core")
    x, y = np.meshgrid(np.linspace(0, 127, 128), np.linspace(0, 127, 128))
    plt.scatter(x, y, c=path_arrival_rates, cmap="YlOrRd", s=0.4)
    plt.colorbar()

    import matplotlib.ticker as ticker
    def create_subplots(c, title=""):
        fig, ax = plt.subplots(nrows=3, ncols=4)
        cmin = 0
        cmax = np.max(c)
        plt.suptitle(title)

        for i in range(0, 12):
            x, y = np.meshgrid(np.linspace(0, 7, 8, dtype=int),
                            np.linspace(0, 3, 4, dtype=int))
            pcm = ax[i//4,i%4].scatter(x, y, c=c[:,:,i], cmap="YlOrRd",
                                       vmin=cmin, vmax=cmax)
            ax[i//4,i%4].set_title(router_link_names[i])
            ax[i//4,i%4].yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

            fig.colorbar(pcm, ax=ax[i//4,i%4])

    #print(router_link_counts[:,:,4])
    #print(router_link_arrival_latency)
    #print(router_link_rates)
    create_subplots(router_link_arrival_rates, "Link Arrival Rates")
    create_subplots(prob_link_blocking, "Probability of Blocking")
    create_subplots(mean_link_waiting_time, "Mean Link Wait Time")
    create_subplots(messages_buffered, "Messaged Buffered")
    create_subplots(router_link_counts, "Link Counts")

    # Plot a heat map
    plt.figure()
    plt.title("Flow Latencies")
    plt.xlabel("Source Core")
    plt.ylabel("Destination Core")
    x, y = np.meshgrid(np.linspace(0, 127, 128), np.linspace(0, 127, 128))
    plt.scatter(x, y, c=flow_latencies, cmap="YlOrRd", s=0.4)
    plt.colorbar()

    plt.figure()
    plt.title("Path Mean Receive Latencies")
    plt.xlabel("Source Core")
    plt.ylabel("Destination Core")
    x, y = np.meshgrid(np.linspace(0, 127, 128), np.linspace(0, 127, 128))
    plt.scatter(x, y, c=path_server_mean_latencies, cmap="YlOrRd", s=0.4)
    plt.colorbar()
    """

    #print(f"mean link transfer delay: {flow_latencies}")
    return flow_latencies, path_counts, path_server_mean_latencies


def schedule_messages_simple(messages):
    neuron_processing = np.zeros((128,))
    message_receiving = np.zeros((128,))
    t = 0.0
    for m in messages:
        neuron_processing[m.src_core] += m.generation_delay
        message_receiving[m.dest_core] += m.receive_delay

    for i in range(0, 128):
        t = max((neuron_processing[i], message_receiving[i], t))
    return t


def schedule_messages_detailed(messages):
    priority = []
    receiving_pipeline_busy = np.zeros((128,))
    # Setup the priority queue, using just a normal list that will be sorted
    #  after every insert. This was easier than figuring out the PriorityQueue
    #  class for floating point priority...
    t = 0.0
    for core in range(0, 128):
        if len(messages[core]) > 0:
            first_message = messages[core][-1]
            priority.append((core, first_message.generation_delay))
    priority.sort(key=lambda x: x[1], reverse=True)

    while len(priority) > 0:
        core, t = priority.pop()
        m = messages[core].pop()
        earliest_receive_time = max(t + m.network_delay,
                                    receiving_pipeline_busy[m.dest_core])
        receiving_pipeline_busy[m.dest_core] = \
            earliest_receive_time + m.receive_delay

        # Push next message
        if len(messages[core]) > 0:
            next_message = messages[core][-1]
            priority.append((core, t+next_message.generation_delay))
            priority.sort(key=lambda x: x[1], reverse=True)

    max_receive = np.max(receiving_pipeline_busy)
    return max(t, max_receive)


class Message:
    """Spike message timing information"""
    def __init__(self, generation_delay=None, receive_delay=None,
                 src_hw=None, dest_hw=None, hops=None, id=None):
        self.id = id
        self.generation_delay = float(generation_delay)
        self.receive_delay = float(receive_delay)
        self.src_tile, self.src_core = hw_str_to_core(src_hw)
        self.dest_tile, self.dest_core = hw_str_to_core(dest_hw)
        self.network_delay = None
        self.hops = hops
        self.pos = None, None, None
        self.pending = False
        self.t_finished_processing = None
        return

    def __repr__(self):
        return str(self)

    def __str__(self):
        s = f"0x{self.id:x}:({self.src_core}->{self.dest_core})@{self.pos}"
        return s


# 1. Read in the network
filename = "dvs_messages.trace"
#filename = "latin_messages.trace"
df = pd.read_csv(filename, converters={"src_hw": str, "dest_hw": str})

timesteps = 128
max_latencies = np.zeros((timesteps,))
mean_latencies = np.zeros((timesteps,))
total_flow_latencies = np.zeros((timesteps,))
total_core_latencies = np.zeros((timesteps,))
max_synapse_processing = np.zeros((timesteps,))
max_neuron_processing = np.zeros((timesteps,))
scheduled_latency = np.zeros((timesteps,))
flow_delays1 = np.zeros((timesteps,))
event_based_latencies = np.zeros((timesteps,))
path_counts = np.zeros((128, 128), dtype=int)  # [src core, dest core]

# Queue of messages for each core
messages = [[] for _ in range(0, 128)]

for timestep in range(0, timesteps):
#for timestep in range(103, 104):
#"""
#for timestep in range(60, 70):
#for timestep in range(62, 63):
    message_generation_latencies = np.zeros((128, 128))
    message_receive_latencies = np.zeros((128, 128))
    df_timestep = df[df["timestep"] == timestep]

    for id, row in df_timestep.iterrows():
        src_tile, src_core = hw_str_to_core(row["src_hw"])
        dest_tile, dest_core = hw_str_to_core(row["dest_hw"])
        messages[src_core].append(Message(row["generation_delay"],
                                       row["processing_latency"],
                                       row["src_hw"], row["dest_hw"],
                                       row["hops"], id))


        message_generation_latencies[src_core, dest_core] += row["generation_delay"]
        """
        if message_generation_latencies[src_core, dest_core] == 0:
            message_generation_latencies[src_core, dest_core] = \
                row["generation_delay"]
        else:
            message_generation_latencies[src_core, dest_core] += \
                min(message_generation_latencies[src_core, dest_core],
                    row["generation_delay"])
        """
        message_receive_latencies[src_core, dest_core] += row["processing_latency"]
        path_counts[src_core, dest_core] += 1

    print(f"Scheduling messages for timestep:{timestep}")
    event_based_latencies[timestep] = sim_schedule_event_based(messages.copy())

    #TODO explore
    #message_generation_latencies = df_timestep["generation_delay"].min()
    #message_generation_latencies = df_timestep["generation_delay"].min()

    flows = np.argwhere(path_counts >= 1)
    flow_delays1[timestep] = np.max(sim_delay_hops(flows))
    flow_delays, flow_counts, _ = sim_delay_mm1k(df_timestep)

    # For all messages, update their network delay based on the average flow
    #  delay
    for core in range(0, 128):
        for m in messages[core]:
            m.network_delay = flow_delays[core, m.dest_core]
    scheduled_latency[timestep] = schedule_messages_detailed(messages)

    print(f"i:{timestep} max flow delay: {np.max(flow_delays):e}")
    max_latencies[timestep] = np.max(flow_delays)

    #print(f"mean server:{message_receive_latencies}")
    #print(f"counts:{flow_counts}")

    max_neuron_processing[timestep] = np.max(np.sum(message_generation_latencies, axis=1))
    max_synapse_processing[timestep] = np.max(np.sum(message_receive_latencies, axis=0))
    max_core = np.argmax(np.sum(message_receive_latencies, axis=0))
    print(f"Max synapse processing happens on core:{max_core}")
    print(message_receive_latencies[:,max_core])
    print(np.sum(message_receive_latencies[:,max_core]))
    #exit()

np.savetxt("runs/analysis/event_based_latencies.csv", event_based_latencies,
           delimiter=",")

#"""

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
LOIHI_TIME_DATA_FILENAME = "loihi_gesture_32x32_time.csv"
NETWORK_DIR = os.path.join(PROJECT_DIR, "runs", "dvs", "loihi_gesture_32x32")
DVS_RUN_DIR = os.path.join(PROJECT_DIR, "runs", "dvs")
LOIHI_TIME_DATA_PATH = os.path.join(DVS_RUN_DIR, LOIHI_TIME_DATA_FILENAME)

np.loadtxt("runs/analysis/event_based_latencies.csv", delimiter=",")
plt.figure()
loihi_data = pd.read_csv(LOIHI_TIME_DATA_PATH)
loihi_times = np.array(loihi_data.loc[:, :] / 1.0e6)
plt.plot(np.arange(1, timesteps-1), loihi_times[0:(timesteps-2), 0] * 1.0e6, "-")
#plt.plot(max_latencies[1:] * 1.0e6)
plt.plot(max_synapse_processing[1:] * 1.0e6, "--")
plt.plot(max_neuron_processing[1:] * 1.0e6, "--")
#plt.plot(scheduled_latency[1:] * 1.0e6)
#plt.plot(flow_delays1[1:] * 1.0e6)
plt.plot(event_based_latencies[1:] * 1.0e6)

#plt.legend(("Measured", "Max", "Mean", "Max Synapse", "Max Neuron"), fontsize=7)
plt.legend(("Measured", "Synapse", "Neuron", "Event-based simulation"))
plt.ylabel("Time-step Latency ($\mu$s)")
plt.xlabel("Time-step")
plt.yticks(np.arange(0, 61, 10))
plt.savefig("runs/dvs_gesture_message_analysis.pdf")

plt.show()
