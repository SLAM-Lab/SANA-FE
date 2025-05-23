// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// fwd.hpp: Forward declarations
#ifndef FWD_HEADER_INCLUDED_
#define FWD_HEADER_INCLUDED_

namespace sanafe
{
class SpikingChip;
struct Timestep;
struct RunData;
struct NocInfo;

struct Scheduler;

class Architecture;
struct NetworkOnChipConfiguration;
struct TileConfiguration;
struct TilePowerMetrics;
struct CoreConfiguration;
struct CorePipelineConfiguration;
struct AxonInConfiguration;
struct PipelineUnitConfiguration;
struct AxonOutConfiguration;
struct AxonOutPowerMetrics;
struct AxonInPowerMetrics;

class SpikingNetwork;
class Neuron;
class NeuronGroup;
struct Connection;
class MappedNeuron;
class MappedConnection;
struct NeuronConfiguration;
struct ModelInfo;
struct PipelineResult;

class Tile;
class Core;
class AxonInUnit;
class AxonOutUnit;
class PipelineUnit;

struct AxonInModel;
struct AxonOutModel;

struct Message;

enum BufferPosition : int;
enum NeuronStatus : uint8_t;

}

// External library forward declarations i.e., not in SANA-FE's namespace
class BookSimConfig;

#endif