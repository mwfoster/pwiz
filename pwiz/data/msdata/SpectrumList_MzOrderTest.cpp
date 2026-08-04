//
// $Id$
//
//
// Original author: Brian Pratt <bspratt .at. proteinms.net>
// AI assistance: Claude Code (Claude Opus 5) <noreply .at. anthropic.com>
//
// Copyright 2026 University of Washington - Seattle, WA
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//


#include "SpectrumList_MzOrder.hpp"
#include "pwiz/utility/misc/unit.hpp"
#include "pwiz/utility/misc/Std.hpp"


using namespace pwiz::cv;
using namespace pwiz::msdata;
using namespace pwiz::util;


namespace {

/// A spectrum with the given peaks, in the given order, plus optionally a third per-peak array
/// standing in for a combined ion mobility axis.
SpectrumPtr makeSpectrum(const string& id, const vector<double>& mzs, const vector<double>& intensities,
                         bool withIonMobilityArray = false)
{
    SpectrumPtr s(new Spectrum);
    s->id = id;
    s->index = 0;
    s->set(MS_MS1_spectrum);
    s->set(MS_ms_level, 1);
    s->setMZIntensityArrays(mzs, intensities, MS_number_of_detector_counts);
    if (withIonMobilityArray)
    {
        BinaryDataArrayPtr ionMobility(new BinaryDataArray);
        ionMobility->set(MS_mean_inverse_reduced_ion_mobility_array);
        vector<double> ionMobilities(mzs.size(), 1.0);
        ionMobility->data.assign(ionMobilities.begin(), ionMobilities.end());
        s->binaryDataArrayPtrs.push_back(ionMobility);
    }
    s->defaultArrayLength = mzs.size();
    return s;
}


SpectrumListSimplePtr makeList(const vector<SpectrumPtr>& spectra)
{
    SpectrumListSimplePtr sl(new SpectrumListSimple);
    for (size_t i = 0; i < spectra.size(); ++i)
    {
        spectra[i]->index = i;
        sl->spectra.push_back(spectra[i]);
    }
    return sl;
}


vector<double> mzsOf(const SpectrumPtr& s)
{
    const BinaryData<double>& data = s->getMZArray()->data;
    return vector<double>(data.begin(), data.end());
}


vector<double> intensitiesOf(const SpectrumPtr& s)
{
    const BinaryData<double>& data = s->getIntensityArray()->data;
    return vector<double>(data.begin(), data.end());
}


// The defect this exists for: peaks stored in ascending intensity rather than ascending m/z.
// Every intensity must still come back attached to the m/z it arrived with - sorting the m/z
// axis alone would leave each value plausible and each pairing wrong, a worse failure than the
// unsorted input.
void testUnsortedSpectrumIsReordered()
{
    vector<double> mzs, intensities;
    mzs.push_back(500.5); intensities.push_back(10);
    mzs.push_back(300.3); intensities.push_back(20);
    mzs.push_back(700.7); intensities.push_back(30);
    mzs.push_back(200.2); intensities.push_back(40);

    vector<SpectrumPtr> spectra;
    spectra.push_back(makeSpectrum("scan=1", mzs, intensities));

    SpectrumList_MzOrder sl(makeList(spectra));
    SpectrumPtr s = sl.spectrum(0, true);

    vector<double> outMzs = mzsOf(s);
    vector<double> outIntensities = intensitiesOf(s);
    unit_assert_operator_equal(4, outMzs.size());
    unit_assert_equal(200.2, outMzs[0], 1e-9); unit_assert_equal(40, outIntensities[0], 1e-9);
    unit_assert_equal(300.3, outMzs[1], 1e-9); unit_assert_equal(20, outIntensities[1], 1e-9);
    unit_assert_equal(500.5, outMzs[2], 1e-9); unit_assert_equal(10, outIntensities[2], 1e-9);
    unit_assert_equal(700.7, outMzs[3], 1e-9); unit_assert_equal(30, outIntensities[3], 1e-9);

    unit_assert(sl.reorderedAnySpectrum());
}


// A file already in m/z order must come through untouched, and must not report itself repaired.
void testSortedSpectrumIsUntouched()
{
    vector<double> mzs, intensities;
    for (int i = 0; i < 20; ++i)
    {
        mzs.push_back(100.0 + i);
        intensities.push_back(20 - i);
    }

    vector<SpectrumPtr> spectra;
    spectra.push_back(makeSpectrum("scan=1", mzs, intensities));

    SpectrumList_MzOrder sl(makeList(spectra));
    SpectrumPtr s = sl.spectrum(0, true);

    unit_assert(mzsOf(s) == mzs);
    unit_assert(intensitiesOf(s) == intensities);
    unit_assert(!sl.reorderedAnySpectrum());
}


// The case a first-spectrum-only probe gets wrong. A short leading spectrum that happens to
// ascend says nothing about the writer - early scans can precede the sample and carry almost no
// peaks - so the checking has to continue until a spectrum with enough peaks settles it.
void testShortOrderedLeaderDoesNotSettleTheFile()
{
    vector<double> shortMzs, shortIntensities;
    shortMzs.push_back(150.1); shortIntensities.push_back(5);
    shortMzs.push_back(250.2); shortIntensities.push_back(7);
    shortMzs.push_back(350.3); shortIntensities.push_back(9);

    vector<double> mzs, intensities;
    mzs.push_back(500.5); intensities.push_back(10);
    mzs.push_back(300.3); intensities.push_back(20);
    mzs.push_back(200.2); intensities.push_back(30);

    vector<SpectrumPtr> spectra;
    spectra.push_back(makeSpectrum("scan=1", shortMzs, shortIntensities));
    spectra.push_back(makeSpectrum("scan=2", mzs, intensities));

    SpectrumList_MzOrder sl(makeList(spectra));
    unit_assert(mzsOf(sl.spectrum(0, true)) == shortMzs);

    vector<double> outMzs = mzsOf(sl.spectrum(1, true));
    unit_assert_equal(200.2, outMzs[0], 1e-9);
    unit_assert_equal(300.3, outMzs[1], 1e-9);
    unit_assert_equal(500.5, outMzs[2], 1e-9);
}


// Once a file is condemned it stays condemned, so a later spectrum that happens to ascend cannot
// switch the checking off and let the ones after it through in writer order.
void testCondemnedFileStaysCondemned()
{
    vector<double> unsortedMzs, unsortedIntensities;
    unsortedMzs.push_back(300.3); unsortedIntensities.push_back(10);
    unsortedMzs.push_back(100.1); unsortedIntensities.push_back(20);

    vector<double> longSortedMzs, longSortedIntensities;
    for (int i = 0; i < 20; ++i)
    {
        longSortedMzs.push_back(100.0 + i);
        longSortedIntensities.push_back(1);
    }

    vector<SpectrumPtr> spectra;
    spectra.push_back(makeSpectrum("scan=1", unsortedMzs, unsortedIntensities));
    spectra.push_back(makeSpectrum("scan=2", longSortedMzs, longSortedIntensities));
    spectra.push_back(makeSpectrum("scan=3", unsortedMzs, unsortedIntensities));

    SpectrumList_MzOrder sl(makeList(spectra));
    sl.spectrum(0, true);
    sl.spectrum(1, true);

    vector<double> outMzs = mzsOf(sl.spectrum(2, true));
    unit_assert_equal(100.1, outMzs[0], 1e-9);
    unit_assert_equal(300.3, outMzs[1], 1e-9);
}


// A combined ion mobility spectrum is legitimately ordered by m/z only within each mobility bin,
// so it must be left exactly as it is - a global sort would shred the bin structure, and the
// third array does not travel with the other two in any case.
void testCombinedIonMobilitySpectrumIsLeftAlone()
{
    vector<double> mzs, intensities;
    mzs.push_back(500.5); intensities.push_back(10);
    mzs.push_back(600.6); intensities.push_back(20);
    mzs.push_back(200.2); intensities.push_back(30);  // roll-over into the next mobility bin
    mzs.push_back(300.3); intensities.push_back(40);

    vector<SpectrumPtr> spectra;
    spectra.push_back(makeSpectrum("merged=1", mzs, intensities, true));

    SpectrumList_MzOrder sl(makeList(spectra));
    SpectrumPtr s = sl.spectrum(0, true);

    unit_assert(mzsOf(s) == mzs);
    unit_assert(intensitiesOf(s) == intensities);
    unit_assert(!sl.reorderedAnySpectrum());
}


// The wrapper is applied to every file pwiz reads, so it has to stay invisible to a file that
// needs nothing. Inventing a DataProcessing would be written out as
// spectrumList/@defaultDataProcessingRef and would move every reference file in the repository.
void testNoDataProcessingIsInvented()
{
    vector<double> mzs, intensities;
    mzs.push_back(100.1); intensities.push_back(1);

    vector<SpectrumPtr> spectra;
    spectra.push_back(makeSpectrum("scan=1", mzs, intensities));

    SpectrumListSimplePtr inner = makeList(spectra);
    unit_assert(!inner->dataProcessingPtr().get());

    SpectrumList_MzOrder sl(inner);
    unit_assert(!sl.dataProcessingPtr().get());
}


// Consumers ask the outermost wrapper whether worker threads are worth it (msconvert does), and
// before this wrapper existed a plain file or vendor list was not a wrapper at all, so the
// question went unanswered and threads stayed on. Inserting this one must not change that.
void testWorkerThreadAnswerIsUnchanged()
{
    vector<double> mzs, intensities;
    mzs.push_back(100.1); intensities.push_back(1);

    vector<SpectrumPtr> spectra;
    spectra.push_back(makeSpectrum("scan=1", mzs, intensities));

    SpectrumList_MzOrder sl(makeList(spectra));
    unit_assert(sl.benefitsFromWorkerThreads());
}

} // namespace


int main(int argc, char* argv[])
{
    TEST_PROLOG(argc, argv)

    try
    {
        testUnsortedSpectrumIsReordered();
        testSortedSpectrumIsUntouched();
        testShortOrderedLeaderDoesNotSettleTheFile();
        testCondemnedFileStaysCondemned();
        testCombinedIonMobilitySpectrumIsLeftAlone();
        testNoDataProcessingIsInvented();
        testWorkerThreadAnswerIsUnchanged();
    }
    catch (exception& e)
    {
        TEST_FAILED(e.what())
    }
    catch (...)
    {
        TEST_FAILED("Caught unknown exception.")
    }

    TEST_EPILOG
}
