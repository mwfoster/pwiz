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


#define PWIZ_SOURCE

#include "SpectrumList_MzOrder.hpp"
#include "pwiz/utility/misc/Std.hpp"
#include <algorithm>


namespace pwiz {
namespace msdata {


namespace {

// A spectrum needs more peaks than this before finding it in m/z order is taken as evidence
// about the writer rather than as coincidence.
const size_t TOO_FEW_PEAKS_TO_MEAN_ANYTHING = 10;

} // namespace


PWIZ_API_DECL SpectrumList_MzOrder::SpectrumList_MzOrder(const SpectrumListPtr& inner)
:   SpectrumListWrapper(inner), verdict_(Verdict::unsettled), reordered_(false)
{
    // SpectrumListWrapper invents a DataProcessing when the inner list has none, and its id is
    // written out as spectrumList/@defaultDataProcessingRef. This wrapper goes on every file read,
    // so inventing one would change the output of every file pwiz writes. It transforms nothing
    // unless a writer got the order wrong, so it stays invisible.
    if (!inner->dataProcessingPtr().get())
        setDataProcessingPtr(DataProcessingPtr());
}


PWIZ_API_DECL SpectrumPtr SpectrumList_MzOrder::spectrum(size_t index, bool getBinaryData) const
{
    SpectrumPtr result = inner_->spectrum(index, getBinaryData);
    ensureMzAscending(result);
    return result;
}


PWIZ_API_DECL SpectrumPtr SpectrumList_MzOrder::spectrum(size_t index, DetailLevel detailLevel) const
{
    SpectrumPtr result = inner_->spectrum(index, detailLevel);
    ensureMzAscending(result);
    return result;
}


PWIZ_API_DECL bool SpectrumList_MzOrder::benefitsFromWorkerThreads() const
{
    // Applied on every file read, so this wrapper must not become the thing that answers a
    // question nothing was answering before. A consumer asks the outermost wrapper (msconvert
    // does); before this wrapper existed, a plain vendor or file list was not a wrapper at all,
    // the question went unanswered and worker threads stayed on. Report that same answer, and
    // defer to a real wrapper when there is one underneath.
    SpectrumListWrapperPtr innerAsWrapper = boost::dynamic_pointer_cast<SpectrumListWrapper>(inner_);
    return innerAsWrapper.get() != NULL ? innerAsWrapper->benefitsFromWorkerThreads() : true;
}


void SpectrumList_MzOrder::ensureMzAscending(const SpectrumPtr& spectrum) const
{
    if (!spectrum.get() || verdict_.load() == Verdict::writerSortsByMz)
        return;

    // Only a plain m/z and intensity spectrum can be reordered on the m/z axis alone. One
    // carrying a third per-peak array - combined ion mobility, or a scanning quadrupole position -
    // is blocked by that axis and ascends in m/z only within each block, so a global sort would
    // destroy it rather than repair it. Such a spectrum is not evidence about the writer either.
    if (spectrum->binaryDataArrayPtrs.size() != 2)
        return;

    BinaryDataArrayPtr mzArray = spectrum->getMZArray();
    BinaryDataArrayPtr intensityArray = spectrum->getIntensityArray();
    if (!mzArray.get() || !intensityArray.get())
        return;

    auto& mzs = mzArray->data;
    auto& intensities = intensityArray->data;
    if (mzs.size() != intensities.size() || mzs.size() < 2)
        return;

    if (std::is_sorted(mzs.begin(), mzs.end()))
    {
        // Peaks in order prove nothing on their own, so the writer is only trusted on the
        // strength of a spectrum with enough peaks to mean it.
        Verdict expected = Verdict::unsettled;
        if (mzs.size() > TOO_FEW_PEAKS_TO_MEAN_ANYTHING)
            verdict_.compare_exchange_strong(expected, Verdict::writerSortsByMz);
        return;
    }

    // One spectrum out of order is proof at any peak count, and it always wins: a writer cannot be
    // talked back into good standing by a later spectrum that happens to ascend.
    verdict_.store(Verdict::writerDoesNotSortByMz);

    vector<pair<double, double> > peaks(mzs.size());
    for (size_t i = 0; i < mzs.size(); ++i)
        peaks[i] = make_pair(mzs[i], intensities[i]);
    std::sort(peaks.begin(), peaks.end(),
              [](const pair<double, double>& a, const pair<double, double>& b) {return a.first < b.first;});
    for (size_t i = 0; i < peaks.size(); ++i)
    {
        mzs[i] = peaks[i].first;
        intensities[i] = peaks[i].second;
    }

    if (!reordered_.exchange(true))
        warn_once(("[SpectrumList_MzOrder] peaks were not written in ascending m/z order (first seen at \"" +
                   spectrum->id + "\"); reordering them. Every consumer assumes that order, so a file "
                   "written in another one extracts as empty without reporting anything.").c_str());
}


} // namespace msdata
} // namespace pwiz
