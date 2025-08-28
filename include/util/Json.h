#pragma once

#include "harpy/Global.h"

#include "rapidjson/document.h"

using RJValue = const rapidjson::Value;
using RJDoc = rapidjson::Document;

void read_json( RJDoc & doc, const string & fn );
