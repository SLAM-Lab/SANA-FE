// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// fwd.hpp: Forward declarations
#ifndef FWD_HEADER_INCLUDED_
#define FWD_HEADER_INCLUDED_

#include <cstdint>

namespace sanafe
{
class SpikingChip;
struct Timestep;
class TimestepHandle;
struct RunData;
class NocInfo;

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
struct NeuronAddress;

class Tile;
class Core;
class AxonInUnit;
class AxonOutUnit;
class PipelineUnit;
class SynapseUnit;
class DendriteUnit;
class SomaUnit;


struct AxonInModel;
struct AxonOutModel;

struct Message;
class CompareMessagesBySentTime;
class CompareMessagesByID;

enum BufferPosition : uint8_t;
enum NeuronStatus : uint8_t;

struct ModelAttribute;

template <typename T>
struct LookupTable;
}

// External library forward declarations i.e., not in SANA-FE's namespace
class BookSimConfig;

#endif