#include "EventFilterManagerGroup.h"

EventFilterManagerGroup::EventFilterManagerGroup() :
	objectsMap(),
	filtersMap()
{
	QString currentFilter("");
	QString currentWatchedObject("");
}

EventFilterManagerGroup::~EventFilterManagerGroup()
{
}

void EventFilterManagerGroup::addObjectToWatch(QString name, QObject* object)
{
	EventFilterManager* eventFilterManager = new EventFilterManager();
	object->installEventFilter(eventFilterManager);
	objectsMap[name] = eventFilterManager;
}

void EventFilterManagerGroup::removeObjectFromWatch(QString name)
{
	EventFilterManager* eventFilterManager = objectsMap[name];
	eventFilterManager->clear();
	eventFilterManager->~EventFilterManager();
	objectsMap.remove(name);
	if (currentWatchedObject == name) {
		currentWatchedObject = "";
		currentFilter = "";
	}
	QList<QString>* filtersToRemove = new QList<QString>();
	QMap<QString, QString>::iterator imap;
	for (imap = filtersMap.begin(); imap != filtersMap.end(); ++imap) {
		if (imap.value() == name)
			filtersToRemove->append(imap.key());
	}
	QList<QString>::const_iterator ilist;
	for (ilist = filtersToRemove->cbegin(); ilist != filtersToRemove->cend(); ++ilist)
	{
		filtersMap.remove(ilist.i->t());
	}
}

void EventFilterManagerGroup::addFilterWidget(const QString targetName, const QString name, QObject* filter)
{
	filtersMap[name] = targetName;
	objectsMap[targetName]->addFilterWidget(name, filter);
}

QObject* EventFilterManagerGroup::removeFilterWidget(const QString name)
{
	QObject* filter = objectsMap[filtersMap[name]]->removeFilterWidget(name);
	filtersMap.remove(name);
	return filter;
}

QObject* EventFilterManagerGroup::getFilterWidget(const QString name)
{
	return objectsMap[filtersMap[name]]->getFilterWidget(name);
}

QString const EventFilterManagerGroup::getCurrentFilterName() const
{
	return currentFilter;
}

void EventFilterManagerGroup::setCurrentFilterName(const QString name)
{
	if (filtersMap.contains(name)) {
		QString newObject = filtersMap[name];
		if (newObject != currentWatchedObject)
			objectsMap[currentWatchedObject]->setCurrentFilterName("");
		objectsMap[newObject]->setCurrentFilterName(name);
		currentWatchedObject = newObject;
		currentFilter = name;
	}
	else if (name == "") {
		if (currentWatchedObject != "")
			objectsMap[currentWatchedObject]->setCurrentFilterName("");
		currentWatchedObject = "";
		currentFilter = "";
	}
	else {
		qDebug() << Q_FUNC_INFO << "The name \"" + name + "\" is not a name of any filter and is not an empty name";
	}
}

QList<QString>* EventFilterManagerGroup::getFilterNames()
{
	return &filtersMap.keys();
}

void EventFilterManagerGroup::clear()
{
	QList<QString>::const_iterator i;
	QList<QString> keys = objectsMap.keys();
	for (i = keys.cbegin(); i != keys.cend(); ++i) {
		objectsMap[*i]->clear();
		objectsMap[*i]->~EventFilterManager();
		objectsMap.remove(*i);
	}
	objectsMap.clear();
	filtersMap.clear();
}

