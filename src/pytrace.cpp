// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#include <ios>
#include <map>
#include <pybind11/gil.h>
#include <pybind11/pytypes.h>
#include <sstream>
#include <stdexcept>
#include <string>

#include "message.hpp"
#include "pytrace.hpp"
#include "timestep.hpp"

pybind11::dict message_to_dict(const sanafe::Message &m)
{
    pybind11::dict message_dict;

    message_dict["generation_delay"] = m.generation_delay;
    message_dict["network_delay"] = m.network_delay;
    message_dict["receive_delay"] = m.receive_delay;
    message_dict["blocked_delay"] = m.blocked_delay;
    message_dict["sent_timestamp"] = m.sent_timestamp;
    message_dict["received_timestamp"] = m.received_timestamp;
    message_dict["processed_timestamp"] = m.processed_timestamp;
    message_dict["timestep"] = m.timestep;
    message_dict["mid"] = m.mid;
    message_dict["spikes"] = m.spikes;
    message_dict["hops"] = m.hops;
    message_dict["src_neuron_offset"] = m.src_neuron_offset;
    message_dict["src_neuron_group_id"] = m.src_neuron_group_id;
    message_dict["src_x"] = m.src_x;
    message_dict["dest_x"] = m.dest_x;
    message_dict["src_y"] = m.src_y;
    message_dict["dest_y"] = m.dest_y;
    message_dict["src_tile_id"] = m.src_tile_id;
    message_dict["src_core_id"] = m.src_core_id;
    message_dict["src_core_offset"] = m.src_core_offset;
    message_dict["dest_tile_id"] = m.dest_tile_id;
    message_dict["dest_core_id"] = m.dest_core_id;
    message_dict["dest_core_offset"] = m.dest_core_offset;
    message_dict["dest_axon_hw"] = m.dest_axon_hw;
    message_dict["dest_axon_id"] = m.dest_axon_id;
    message_dict["placeholder"] = m.placeholder;

    return message_dict;
}

std::map<std::string, PerfStatistic> timestep_data_to_map(
        const sanafe::Timestep &ts)
{
    std::map<std::string, PerfStatistic> ts_map;

    ts_map["timestep"] = ts.timestep;
    ts_map["fired"] = ts.neurons_fired;
    ts_map["updated"] = ts.neurons_updated;
    ts_map["hops"] = ts.total_hops;
    ts_map["spikes"] = ts.spike_count;
    ts_map["sim_time"] = ts.sim_time;
    // Energy values
    ts_map["synapse_energy"] = ts.synapse_energy;
    ts_map["dendrite_energy"] = ts.dendrite_energy;
    ts_map["soma_energy"] = ts.soma_energy;
    ts_map["network_energy"] = ts.network_energy;
    ts_map["total_energy"] = ts.total_energy;

    return ts_map;
}

void PyTrace::open(pybind11::object trace_obj, const bool overwrite_trace)
{
    mode = TraceMode::trace_none;
    // The trace object can be None, a boolean (True/False) enabling memory
    //  based traces, a string specifying a trace filename or a Python file
    //  object
    obj = trace_obj;

    if (!trace_obj.is_none())
    {
        // Check if it's a boolean (True)
        if (pybind11::isinstance<pybind11::bool_>(trace_obj))
        {
            const bool trace_enabled = trace_obj.cast<bool>();
            if (trace_enabled)
            {
                mode = TraceMode::trace_memory;
            }
        }
        // Check if it's a file-like object (i.e., has a 'write' method)
        else if (pybind11::hasattr(trace_obj, "write"))
        {
            mode = TraceMode::trace_file;
            // File handle will be used directly via Python calls

            if (overwrite_trace)
            {
                obj.attr("seek")(0);
            }
        }
        // Check if it's a string (filename)
        else if (pybind11::isinstance<pybind11::str>(trace_obj))
        {
            mode = TraceMode::trace_file;
            const auto filename = trace_obj.cast<std::string>();

            const std::ios::openmode flags = overwrite_trace ?
                    (std::ios::out | std::ios::trunc) :
                    (std::ios::out | std::ios::app);
            file.open(filename, flags);
            if (!file.is_open())
            {
                throw std::runtime_error(
                        "Failed to open trace file: " + filename);
            }
        }
        else
        {
            throw std::invalid_argument(
                    "trace_obj must be None, True, a filename string, "
                    "or a file-like object");
        }
    }
}

void PyTrace::write_header()
{
    if (mode == TraceMode::trace_file)
    {
        std::ostringstream trace_ss;
        get_header_string(trace_ss);
        if (file.is_open())
        {
            // Write to C++ file stream
            file << trace_ss.str();
            file.flush(); // Ensure data is written
        }
        else if (!obj.is_none())
        {
            // Write to Python file handle
            const pybind11::gil_scoped_acquire acquire_for_write;
            obj.attr("write")(trace_ss.str());
            if (pybind11::hasattr(obj, "flush"))
            {
                obj.attr("flush")();
            }
        }
    }
}

void PyHardwareTrace::record_hw_metrics(const sanafe::Timestep &ts)
{
    if (mode == TraceMode::trace_memory)
    {
        // Note we must always acquire the GIL before modifying Python objects
        const pybind11::gil_scoped_acquire acquire;
        // Make call to appropriate overriden virtual routine
        append_hw_metrics_to_memory(ts);
    }
    else if (mode == TraceMode::trace_file)
    {
        std::ostringstream trace_ss;
        // Make call to appropriate overriden virtual routine
        write_hw_trace_to_stream(trace_ss, ts);
        if (file.is_open())
        {
            // Write to C++ file stream
            file << trace_ss.str();
            file.flush(); // Ensure data is written
        }
        else if (!obj.is_none())
        {
            // Write to Python file handle
            const pybind11::gil_scoped_acquire acquire_for_write;
            obj.attr("write")(trace_ss.str());
            if (pybind11::hasattr(obj, "flush"))
            {
                obj.attr("flush")();
            }
        }
    }
    // else do nothing
}

void PyNetworkTrace::record_net_activity(const long int timestep)
{
    if (mode == TraceMode::trace_memory)
    {
        // Note we must always acquire the GIL before modifying Python objects
        const pybind11::gil_scoped_acquire acquire;
        // Make call to appropriate overriden virtual routine
        append_net_activity_to_memory();
    }
    else if (mode == TraceMode::trace_file)
    {
        std::ostringstream trace_ss;
        // Make call to appropriate overriden virtual routine
        write_net_trace_to_stream(trace_ss, timestep);
        if (file.is_open())
        {
            // Write to C++ file stream
            file << trace_ss.str();
            file.flush(); // Ensure data is written
        }
        else if (!obj.is_none())
        {
            // Write to Python file handle
            const pybind11::gil_scoped_acquire acquire_for_write;
            obj.attr("write")(trace_ss.str());
            if (pybind11::hasattr(obj, "flush"))
            {
                obj.attr("flush")();
            }
        }
    }
    // else do nothing
}
