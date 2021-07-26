/*
    File: monomerCoordinates.cc
*/
/*
Open Source License
Copyright (c) 2016, Christian E. Schafmeister
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 
This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
#define	DEBUG_LEVEL_NONE


#include <cando/chem/monomerCoordinates.h>
#include <cando/chem/coordSys.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/superposableConformationCollection.h>
#include <cando/chem/topology.h>
#include <cando/chem/monomer.h>
#include <clasp/core/environment.h>
#include <cando/chem/constitution.h>
#include <clasp/core/lispStream.h>
#include <cando/chem/rotamer.h>
#include <clasp/core/wrappers.h>


namespace chem {




//
// Constructor
//

//
// Destructor
//

void	MonomerCoordinates_O::initialize()
{
    this->_TrainerName = "";
}


#ifdef XML_ARCHIVE
    void	MonomerCoordinates_O::archiveBase(core::ArchiveP node)
{
    node->attribute("comment",this->_Comment);
    node->attributeIfNotDefault<string>("trainerName",this->_TrainerName, "" );
    node->attribute("context",this->_Context);
    node->archiveVector0("rotamers",this->_Rotamers);
}
#endif


  MonomerCoordinates_sp MonomerCoordinates_O::make(MonomerContext_sp context)
  {
    auto  me  = gctools::GC<MonomerCoordinates_O>::allocate_with_default_constructor();
    me->setContext(context);
    return me;
  };






void	MonomerCoordinates_O::_defineFromConformationExplorerOrDebug(
				Topology_sp	topology,
				ConformationExplorer_sp conformationExplorer,
				Monomer_sp focusMonomer,
				bool debug,
				uint maximumNumberOfConformations )
{
    IMPLEMENT_ME(); // update
#if 0
    ExtractFragmentBase_sp	extractFragment;
    SuperposableConformationCollection_sp uniqueStructures;
    FragmentCoordinates_sp	fragmentCoordinates;
    FragmentCoordinates_sp	coreFragmentCoordinates;
    ScaffoldList_sp		scaffoldList;
    ExtractScaffold_sp extractScaffold;
    Aggregate_sp		agg;
    ExtractInternalFrameBase_sp	firstExtractFrame;
    Frame_sp			firstFrame;

    //
    // Extract the core FragmentCoordinates
    // and the ScaffoldList using
    // the same conformationExplorer
    //

    { _BLOCK_TRACE("Extracting ScaffoldList and core FragmentCoordinates");
	{ _BLOCK_TRACE("Identify unique core conformations");
	    extractFragment = topology->getExtractCoreFragment();
	    LOG(BF("Got core fragment: %s") % extractFragment->description().c_str() );
	    uniqueStructures = extractFragment->isolateUniqueStructures(
		conformationExplorer, focusMonomer);

	    extractScaffold = topology->getExtractScaffold();
	    Anchor_sp anchor = extractScaffold->getAnchor();
	    //
	    // If the Anchor is AnchorOrigin then we want to move
	    // the molecule to the origin so it isnt out in space somewhere
	    //
	    if ( anchor.isA<AnchorOrigin_O>() )
	    {
		// This is anchored on the origin so we want to translate
		// every structure in uniqueStructures to the origin
		//
		agg = uniqueStructures->getMatter().as<Aggregate_O>();
		// We are going to transform every structure to the
		// first complete frame in the extract scaffold
		// there has to be at least one.
		// The only time there wouldn't be is if the anchor is
		// not an AnchorOrigin, which it has to be according
		// to the test above
		ASSERT(extractScaffold->numberOfExtractCompleteFrames()>0);
		firstExtractFrame = *(extractScaffold->begin_ExtractCompleteFrames());
		firstFrame = firstExtractFrame->getFrame();
		ConformationCollection_O::entryIterator sli;
		CoordinateSystem_sp	coordSys;
		Matrix			toOriginTransform;
		for ( sli=uniqueStructures->begin_Entries();
		      sli!=uniqueStructures->end_Entries(); sli++ )
		{
		    (*sli)->writeCoordinatesToMatter(agg);
		    TRY()
		    {
		        coordSys = firstFrame->getCoordinateSystem(focusMonomer);
		    } catch (core::Condition& e)
		      {
			  if ( !e.conditionObject().isA<core::MathException_O>() )
			  {
			      RETHROW(e.conditionObject());
			  }
			  if ( !debug )
			  {
			      RETHROW(e.conditionObject());
			  }
			  coordSys = CoordinateSystem_O::create();
		      }
		    toOriginTransform = coordSys->matrixToCanonical();
		    agg->applyTransformToAtoms(toOriginTransform);
		    (*sli)->extractCoordinatesFromMatter(agg);
		}
	    }
	    stringstream problems;
	    if ( !extractFragment->sharesAtomsWithExtractScaffold(extractScaffold,
								  focusMonomer, problems ) )
	    {
		stringstream ss;
		ss << "ExtractFragment does not share atoms with ExtractScaffold"<<std::endl;
		ss << "ExtractFragment constitution["<< topology->getConstitution()->getName() << "] "
		   << "topology[" << topology->getName() << "]"
		   << " fragmentName:("+extractFragment->getFragment()->getName()+") atoms: ["
		   << extractFragment->getFragment()->getAtomNames()->asString() << "]" << std::endl;
		ss << problems.str();
		SIMPLE_ERROR(BF(ss.str()));
	    }
	}
	{ _BLOCK_TRACE("Extract core FragmentCoordinates");
	    coreFragmentCoordinates = FragmentCoordinates_O::create();
	    coreFragmentCoordinates->defineFromConformationCollection(extractFragment,
								      focusMonomer,
								      uniqueStructures, debug, 
								      maximumNumberOfConformations );
	    this->_FragmentCoordinates.clear();
	    this->_FragmentCoordinates.push_back(coreFragmentCoordinates);
	}
	{ _BLOCK_TRACE("Extract ScaffoldList");
	    scaffoldList = ScaffoldList_O::create();
	    scaffoldList->defineFromExtractScaffoldAndConformationCollection( extractScaffold,
									      focusMonomer,
									      uniqueStructures,
									      coreFragmentCoordinates,
									      debug,
									      maximumNumberOfConformations);
	    this->_ScaffoldList = scaffoldList;
	}
    }
    //
    // Extract the other fragments
    //
    {_BLOCK_TRACE("Extract the non-core fragments");
	Topology_O::iterateExtractFragments	fi;
	for ( fi = topology->begin_ExtractFragments();
	      fi!=topology->end_ExtractFragments(); fi++ )
	{
	    uniqueStructures = (*fi)->isolateUniqueStructures(conformationExplorer,
							      focusMonomer);
	    fragmentCoordinates = FragmentCoordinates_O::create();
	    fragmentCoordinates->defineFromConformationCollection(*fi,
								  focusMonomer,
								  uniqueStructures,
								  debug,
								  maximumNumberOfConformations );
	    this->_FragmentCoordinates.push_back(fragmentCoordinates);
	}
    }
    LOG(BF("MonomerCoordinates defined: %s") % this->asXmlString().c_str()  );
#endif
}


CL_LISPIFY_NAME("defineFromConformationExplorer");
CL_DEFMETHOD void	MonomerCoordinates_O::defineFromConformationExplorer(
				Topology_sp	topology,
				ConformationExplorer_sp conformationExplorer,
				Monomer_sp focusMonomer,
				uint maximumNumberOfConformations )
{
    this->_defineFromConformationExplorerOrDebug(topology, conformationExplorer, focusMonomer, false,
    						maximumNumberOfConformations);
}



/*
 * Carry out a MonomerCoordinate extraction test and return any problems that were encountered
 * in a Cons containing 
 * :sawProblems xxx  and :problems yyy pairs
 */

CL_LISPIFY_NAME("testExtraction");
CL_DEFMETHOD core::List_sp	MonomerCoordinates_O::testExtraction(
				Topology_sp	topology,
				Aggregate_sp	aggregate,
				Monomer_sp focusMonomer )
{
    ConformationExplorer_sp conformationExplorer = ConformationExplorer_O::create();
    conformationExplorer->setMatter(aggregate);
    conformationExplorer->superposeAllHeavyAtoms();
    ConformationExplorerEntry_sp entry = conformationExplorer->createEntry();
    ConformationExplorerEntryStage_sp stage = entry->createEntryStage(core::SimpleBaseString_O::make("test"));
    LOG(BF("Setting stage") );
    stage->setExternalInterfaceName(core::SimpleBaseString_O::make("testPackage"));
    stage->setModel(core::SimpleBaseString_O::make("testModel"));
    stage->setEnergyKCal(987654321.0);
    entry->setSelectedStage(stage);
    stage->extractCoordinatesFromMatter(aggregate);
    stage->setComplete(true);
    conformationExplorer->appendEntry(entry);
    core::T_sp problems = core::clasp_make_string_output_stream();
    IMPLEMENT_MEF("The following code catches an exception - figure out how to do it using Conditions");
#if 0
    TRY()
    {
        this->_defineFromConformationExplorerOrDebug(topology,conformationExplorer,focusMonomer,true,UndefinedUnsignedInt);
    } catch ( core::Condition& e )
    {
	if ( !e.conditionObject().isA<core::MathException_O>() )
	{
	    throw;
//	    RESIMPLE_ERROR(BF("An error was caught in MonomerCoordinates_O::testExtraction - here it is: %s") % e.message() ));
	}
	sawProblem = true;
	LOG(BF("%s") % e.message().c_str()  );
	problems->write(e.message());
    }
    core::Cons_sp results = core::Cons_O::createList(_lisp->internKeyword("sawProblems"),
					 _lisp->_boolean(sawProblem),
					 _lisp->internKeyword("problems"),
					 problems);
    return results;
#endif
}














};
